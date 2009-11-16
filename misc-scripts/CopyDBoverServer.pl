#!/usr/local/bin/perl -w

use strict;
use warnings;

use DBI;
use File::Copy;
use File::Spec::Functions
  qw( rel2abs curdir canonpath updir catdir catfile );
use Getopt::Long;
use IO::File;
use Sys::Hostname;

my $start_time = time();

sub short_usage {
  my $indent = ' ' x length($0);

  print <<EOT;
Usage:
  $0 --pass=XXX [--noflush] [--nocheck] [--force] \\
  $indent [--subset=XXX] [--help] input_file\

  Use --help to get a much longer help text.

EOT
}

sub long_usage {
  my $indent = ' ' x length($0);

  print <<EOT;
Usage:
  $0 --pass=XXX [--noflush] [--nocheck] [--force] \\
  $indent [--subset=XXX] [--help] input_file

Description:

  Safetly copy MySQL databases between different servers and run
  myisamchk on the indexes when done.

  The script A) transfers the database files to a local staging
  directory, B) checks the files for errors using myisamchk, and C)
  moves them into place in the database server directory.


Command line switches:

  --pass=XXX        (Required)
                    The password for the 'ensadmin' MySQL user to
                    connect to the database.

  --noflush         (Optional)
                    Skips table flushing completely.  Use very
                    sparsingly as copying un-flushed databases
                    potentially means missing data not yet flushed to
                    disk.  Use only after due consideration.

  --nocheck         (Optional)
                    Skip running myisamchk on the copied table files.
                    Use this only if you are *absolutely* sure that the
                    table files are ok and that you really do not have
                    time to wait for the check to run.  Use only after
                    due consideration.

  --force           (Optional)
                    Ordinarily, the script refuses to overwrite an
                    already existing staging directory.  This switch
                    will force it to re-use the directory if it exists.
                    This is useful for continuing an aborted copy.

                    NOTE: This switch will not ever force overwriting of
                    a server database directory.

  --subset=XXX      (Optional)
                    Only copy a subset of the tables.  Currently, the
                    only subsets defined are 'xref' and 'stable_id'.

  --help            (Optional)
                    Displays this text.


Input file format:

  The input file should be a tab/space-separated file with six fields
  per line:

    1. Source server
    2. Source server port
    3. Source database name

    4. Target server
    5. Target server port
    6. Target database name

  For example:

  genebuild1 3306 at6_gor2_wga ens-staging1 3306 gorilla_gorilla_core_57_1

  (with each field being separated by a tab or space).

  Blank lines, lines containing only whitespaces, and lines starting
  with '#', are silently ignored.


Script restrictions:

  1. You must run the script on the destination server.

  2. The script will only allow copying databases between a preset set
     of servers.

  3. The script must be run as the 'mysqlens' Unix user.  Talk to a
     recent release coordinator for access.

  4. The script will only copy MYISAM tables.  Databases with InnoDB
     tables will have to be copied manually using mysqldump.  InnoDB
     tables will make the script throw an error in the table checking
     stage.

EOT
} ## end sub long_usage

my (
  $opt_password, $opt_flush,  $opt_check,
  $opt_force,    $opt_subset, $opt_help
);

$opt_flush = 1;    # Flush by default.
$opt_check = 1;    # Check tables by default.
$opt_force = 0;    # Do not reuse existing staging directory by default.

if (
  !GetOptions(
    'pass=s'   => \$opt_password,
    'flush!'   => \$opt_flush,
    'check!'   => \$opt_check,
    'force!'   => \$opt_force,
    'subset=s' => \$opt_subset,
    'help'     => \$opt_help
  )
  || ( !defined($opt_password) && !defined($opt_help) ) )
{
  short_usage();
  exit 1;
}

if ( defined($opt_help) ) {
  long_usage();
  exit 0;
}

if ( scalar( getpwuid($<) ) ne 'mysqlens' ) {
  die("You need to run this script as the 'mysqlens' user.\n");
}

my $input_file = shift(@ARGV);

if ( !defined($input_file) ) {
  short_usage();
  exit 1;
}

my %executables = (
  'myisamchk' => '/software/farm/mysql/bin/myisamchk',
  'rsync'     => '/usr/bin/rsync'
);

my %table_subsets = (
  'xref' => [
    'xref',            'object_xref',
    'identity_xref',   'go_xref',
    'external_db',     'external_synonym',
    'unmapped_reason', 'unmapped_object',
  ],
  'stable_id' => [
    'gene_stable_id',        'transcript_stable_id',
    'translation_stable_id', 'exon_stable_id',
    'mapping_session',       'stable_id_event',
    'gene_archive',          'peptide_archive',
  ] );

if ( defined($opt_subset) ) {
  if ( !exists( $table_subsets{ lc($opt_subset) } ) ) {
    die( sprintf( "Unknown table subset '%s'\n", $opt_subset ) );
  }
  $opt_subset = lc($opt_subset);
}

my $run_hostname = ( gethostbyname( hostname() ) )[0];
my $working_dir  = rel2abs( curdir() );

$run_hostname =~ s/\..+//;    # Cut off everything but the first part.

##====================================================================##
##  Read the configuration file line by line and try to validate all  ##
##  parts of each line.  Store the validated in the @todo list (a     ##
##  list of hashes) for later processing.                             ##
##====================================================================##

my $in = IO::File->new( '<' . $input_file )
  or die( sprintf( "Can not open '%s' for reading", $input_file ) );

my @todo;                     # List of verified databases etc. to copy.

my $lineno = 0;
while ( my $line = $in->getline() ) {
  ++$lineno;

  $line =~ s/^\s+//;          # Strip leading whitespace.
  $line =~ s/\s+$//;          # Strip trailing whitespace.

  if ( $line =~ /^\#/ )   { next }    # Comment line.
  if ( $line =~ /^\s*$/ ) { next }    # Empty line.

  my $failed = 0;                     # Haven't failed so far...

  my (
    $source_server, $source_port, $source_db,
    $target_server, $target_port, $target_db
  ) = split( /\s+/, $line );

  my $source_hostname = ( gethostbyname($source_server) )[0];
  my $target_hostname = ( gethostbyname($target_server) )[0];

  # Verify source server and port.
  if ( !defined($source_hostname) || $source_hostname eq '' ) {
    warn(
      sprintf(
        "line %d: Source server '%s' is not valid.\n",
        $lineno, $source_server
      ) );
    $failed = 1;
  } else {
    $source_hostname =~ s/\..+//;
  }

  if ( !defined($source_port) || $source_port =~ /\D/ ) {
    warn(
      sprintf(
        "line %d: Source port '%s' is not a number.\n",
        $lineno, $source_port || ''
      ) );
    $failed = 1;
  }

  # Verify target server and port.
  if ( !defined($target_hostname) || $target_hostname eq '' ) {
    warn(
      sprintf(
        "line %d: Target server '%s' is not valid.\n",
        $lineno, $target_server
      ) );
    $failed = 1;
  } else {
    $target_hostname =~ s/\..+//;
  }

  if ( !defined($target_port) || $target_port =~ /\D/ ) {
    warn(
      sprintf(
        "line %d: Target port '%s' is not a number.\n",
        $lineno, $target_port || ''
      ) );
    $failed = 1;
  }

  # Make sure we running on the target server.
  if ( !$failed && $run_hostname ne $target_hostname ) {
    warn(
      sprintf(
        "line %d: "
          . "This script needs to be run on the destination server "
          . "'%s' ('%s').\n",
        $lineno, $target_server, $target_hostname
      ) );
    $failed = 1;
  }

  if ( !$failed ) {
    push(
      @todo,
      {
        'source_server'   => $source_server,
        'source_hostname' => $source_hostname,
        'source_port'     => $source_port,
        'source_db'       => $source_db,
        'target_server'   => $target_server,
        'target_hostname' => $target_hostname,
        'target_port'     => $target_port,
        'target_db'       => $target_db,
      } );
  }
} ## end while ( my $line = $in->getline...)

$in->close();

##====================================================================##
##  Take the copy specifications from the @todo list and for each     ##
##  specification copy the database to a staging area using rsync,    ##
##  check it with myisamchk, and move it in place in the database     ##
##  directory.                                                        ##
##====================================================================##

foreach my $spec (@todo) {
  my $source_server   = $spec->{'source_server'};
  my $source_hostname = $spec->{'source_hostname'};
  my $source_port     = $spec->{'source_port'};
  my $source_db       = $spec->{'source_db'};
  my $target_server   = $spec->{'target_server'};
  my $target_hostname = $spec->{'target_hostname'};
  my $target_port     = $spec->{'target_port'};
  my $target_db       = $spec->{'target_db'};

  my $label = sprintf( "{ %s -> %s }==", $source_db, $target_db );
  print( '=' x ( 80 - length($label) ), $label, "\n" );

  print("CONNECTING TO SOURCE AND TARGET DATABASES "
      . "TO GET 'datadir'\n" );

  my $source_dsn = sprintf( "DBI:mysql:database=%s;host=%s;port=%d",
    $source_db, $source_hostname, $source_port );

  my $source_dbh =
    DBI->connect( $source_dsn, 'ensadmin', $opt_password,
    { 'PrintError' => 1, 'AutoCommit' => 0 } );

  if ( !defined($source_dbh) ) {
    warn(
      sprintf(
        "Failed to connect to the source database '%s:%d/%s'.\n",
        $source_server, $source_port, $source_db
      ) );

    $spec->{'status'} =
      sprintf( "FAILED: can not connect to source database '%s:%d/%s'.",
      $source_server, $source_port, $source_db );
    next;
  }

  my $target_dsn = sprintf( "DBI:mysql:host=%s;port=%d",
    $target_hostname, $target_port );

  my $target_dbh =
    DBI->connect( $target_dsn, 'ensadmin', $opt_password,
    { 'PrintError' => 1, 'AutoCommit' => 0 } );

  if ( !defined($target_dbh) ) {
    warn(
      sprintf(
        "Failed to connect to the target server '%s:%d'.\n",
        $target_server, $target_port
      ) );

    $spec->{'status'} =
      sprintf( "FAILED: can not connect to target server '%s:%d'.",
      $target_server, $target_port );

    $source_dbh->disconnect();
    next;
  }

  # Get source and target server data directories.
  my $source_dir =
    $source_dbh->selectall_arrayref("SHOW VARIABLES LIKE 'datadir'")
    ->[0][1];
  my $target_dir =
    $target_dbh->selectall_arrayref("SHOW VARIABLES LIKE 'datadir'")
    ->[0][1];

  $target_dbh->disconnect();

  if ( !defined($source_dir) ) {
    warn(
      sprintf(
        "Failed to find data directory for source server at '%s:%d'.\n",
        $source_server, $source_port
      ) );

    $spec->{'status'} = sprintf(
      "FAILED: can not find data directory on source server '%s:%d'.",
      $source_server, $source_port );

    $source_dbh->disconnect();
    next;
  }

  if ( !defined($target_dir) ) {
    warn(
      sprintf(
        "Failed to find data directory for target server at '%s:%d'.\n",
        $target_server, $target_port
      ) );

    $spec->{'status'} = sprintf(
      "FAILED: can not find data directory on target server '%s:%d'.",
      $target_server, $target_port );

    $source_dbh->disconnect();
    next;
  }

  printf( "SOURCE 'datadir' = '%s'\n", $source_dir );
  printf( "TARGET 'datadir' = '%s'\n", $target_dir );

  my $tmp_dir = canonpath( catdir( $target_dir, updir(), 'tmp' ) );

  my $staging_dir = catdir( $tmp_dir, sprintf( "tmp.%s", $target_db ) );
  my $destination_dir = catdir( $target_dir, $target_db );

  $spec->{'status'} = 'SUCCESS';    # Assume success until failure.

  # Try to make sure the temporary directory and the final destination
  # directory actually exists, and that the staging directory within the
  # temporary directory does *not* exist.  Allow the staging directory
  # to be reused when the --force switch is used.

  if ( !-d $tmp_dir ) {
    die(
      sprintf( "Can not find the temporary directory '%s'", $tmp_dir )
    );
  }

  if ( -d $destination_dir ) {
    warn(
      sprintf( "Destination directory '%s' already exists.\n",
        $destination_dir ) );

    $spec->{'status'} =
      sprintf( "FAILED: database destination directory '%s' exists.",
      $destination_dir );

    $source_dbh->disconnect();
    next;
  }

  if ( !$opt_force && -d $staging_dir ) {
    warn(
      sprintf( "Staging directory '%s' already exists.\n",
        $staging_dir ) );

    $spec->{'status'} =
      sprintf( "FAILED: staging directory '%s' exists.", $staging_dir );

    $source_dbh->disconnect();
    next;
  }

  if ( !mkdir($staging_dir) ) {
    if ( !$opt_force || !-d $staging_dir ) {
      warn(
        sprintf( "Failed to create staging directory '%s'.\n",
          $staging_dir ) );

      $spec->{'status'} =
        sprintf( "FAILED: can not create staging directory '%s'.",
        $staging_dir );

      $source_dbh->disconnect();
      next;
    }
  }

  my @tables;
  if ( defined($opt_subset) ) {
    # Copy only a subset of the tables.
    @tables = @{ $table_subsets{$opt_subset} };
  } else {
    # Copy all tables.

    foreach
      my $table ( @{ $source_dbh->selectall_arrayref('SHOW TABLES') } )
    {
      push( @tables, $table->[0] );
    }
  }

  # Lock tables with a read lock.
  print("LOCKING TABLES...\n");
  $source_dbh->do(
    sprintf( "LOCK TABLES %s READ", join( ' READ, ', @tables ) ) );

  if ($opt_flush) {
    # Flush tables.

    print("FLUSHING TABLES...\n");
    $source_dbh->do(
      sprintf( "FLUSH TABLES %s", join( ', ', @tables ) ) );
  }

  ##------------------------------------------------------------------##
  ## COPY                                                             ##
  ##------------------------------------------------------------------##

  print( '-' x 37, ' COPY ', '-' x 37, "\n" );

  # Set up database copying.  We're using rsync for this because it's
  # using SSH for network transfers, because it may be used for local
  # copy too, and because it has good inclusion/exclusion filter
  # options.

  my @copy_cmd = ( $executables{'rsync'}, '--archive', '--progress' );

  if ($opt_force) {
    push( @copy_cmd, '--delete', '--delete-excluded' );
  }

  if ( defined($opt_subset) ) {
    # Subset copy.
    push( @copy_cmd, "--include=db.opt" );
    push( @copy_cmd,
      map { sprintf( "--include=%s.*", $_ ) }
        @{ $table_subsets{$opt_subset} } );
    push( @copy_cmd, "--exclude=*" );
  }

  if ( $source_hostname eq $target_hostname ) {
    # Local copy.
    push( @copy_cmd,
      sprintf( "%s/", catdir( $source_dir, $source_db ) ) );
  } else {
    # Copy from remote server.
    push(
      @copy_cmd,
      sprintf( "%s:%s/",
        $source_hostname, catdir( $source_dir, $source_db ) ) );
  }

  push( @copy_cmd, sprintf( "%s/", $staging_dir ) );

  # Perform the copy and make sure it succeeds.

  printf( "COPYING '%s:%d/%s' TO STAGING DIRECTORY '%s'\n",
    $source_server, $source_port, $source_db, $staging_dir );

  # For debugging:
  # print( join( ' ', @copy_cmd ), "\n" );

  my $copy_failed = 0;
  if ( system(@copy_cmd) != 0 ) {
    warn(
      sprintf(
        "Failed to copy database.\n"
          . "Please clean up '%s' (if needed).",
        $staging_dir
      ) );
    $copy_failed = 1;
  }

  # Unlock tables.
  print("UNLOCKING TABLES...\n");
  $source_dbh->do('UNLOCK TABLES');

  $source_dbh->disconnect();

  if ($copy_failed) {
    $spec->{'status'} =
      sprintf( "FAILED: copy failed (cleanup of '%s' may be needed).",
      $staging_dir );
    next;
  }

  ##------------------------------------------------------------------##
  ## CHECK                                                            ##
  ##------------------------------------------------------------------##

  print( '-' x 36, ' CHECK ', '-' x 37, "\n" );

  # Check the copied table files with myisamchk.  Let myisamchk
  # automatically repair any broken or un-closed tables.

  if ( !$opt_check ) {
    print("NOT CHECKING...\n");
  } else {
    print("CHECKING TABLES...\n");

    my $check_failed = 0;

    foreach my $table (@tables) {
      my @check_cmd = (
        $executables{'myisamchk'},
        '--force',
        '--check',
        '--check-only-changed',
        '--update-state',
        '--silent',
        '--silent',    # Yes, twice.
        catfile( $staging_dir, $table ) );

      if ( system(@check_cmd) != 0 ) {
        $check_failed = 1;

        warn(
          sprintf(
            "Failed to check some tables. "
              . "Is this an InnoDB database maybe?\n"
              . "Please clean up '%s'.\n",
            $staging_dir
          ) );

        $spec->{'status'} = sprintf(
          "FAILED: MYISAM table check failed "
            . "(cleanup of '%s' may be needed).",
          $staging_dir
        );

        last;
      }
    } ## end foreach my $table (@tables)

    if ($check_failed) { next }

  } ## end else [ if ( !$opt_check ) ]

  ##------------------------------------------------------------------##
  ## MOVE                                                             ##
  ##------------------------------------------------------------------##

  print( '-' x 37, ' MOVE ', '-' x 37, "\n" );

  # Move table files into place in $destination_dir using
  # File::Copy::move(), and remove the staging directory.  We already
  # know that the destination directory does not exist.

  printf( "MOVING '%s' TO '%s'...\n", $staging_dir, $destination_dir );

  if ( !mkdir($destination_dir) ) {
    warn(
      sprintf(
        "Failed to create destination directory '%s'.\n"
          . "Please clean up '%s'.\n",
        $destination_dir, $staging_dir
      ) );

    $spec->{'status'} =
      sprintf( "FAILED: can not create destination directory '%s' "
        . "(cleanup of '%s' may be needed)",
      $destination_dir, $staging_dir );
    next;
  }

  move( catfile( $staging_dir, 'db.opt' ), $destination_dir );

  foreach my $table ( @tables ) {
    my @files =
      glob( catfile( $staging_dir, sprintf( "%s.*", $table ) ) );

    printf( "Moving %s...\n", $table );

    foreach my $file (@files) {
      if ( !move( $file, $destination_dir ) ) {
        warn(
          sprintf(
            "Failed to move database.\n"
              . "Please clean up '%s' and '%s'.\n",
            $staging_dir, $destination_dir
          ) );

        $spec->{'status'} =
          sprintf( "FAILED: move from staging directory failed "
            . "(cleanup of '%s' and '%s' may be needed)",
          $staging_dir, $destination_dir );
        next;
      }
    }
  }

  # Remove the now empty staging directory.
  if ( !rmdir($staging_dir) ) {
    warn(
      sprintf(
        "Failed to unlink the staging directory '%s'.\n"
          . "Clean this up manually.\n",
        $staging_dir
      ) );

    $spec->{'status'} =
      sprintf( "SUCCESS: cleanup of '%s' may be needed", $staging_dir );
  }

} ## end foreach my $spec (@todo)

# Display summary.

my $label = '{ SUMMARY }~~';
print( "\n", '~' x ( 80 - length($label) ), $label, "\n" );

foreach my $spec (@todo) {
  printf( "%s -> %s\n  %s\n\n",
    $spec->{'source_db'}, $spec->{'target_db'}, $spec->{'status'} );
}

print("DONE!\n\n");

END {
  my $seconds = time() - $start_time;

  my $hours = int( $seconds/( 60*60 ) );
  $seconds -= 60*60*$hours;

  my $minutes = int( $seconds/60 );
  $seconds -= 60*$minutes;

  printf(
    "Time taken: %s%dm%ds\n",
    ( $hours > 0 ? sprintf( "%dh", $hours ) : '' ),
    $minutes, $seconds
  );
}
