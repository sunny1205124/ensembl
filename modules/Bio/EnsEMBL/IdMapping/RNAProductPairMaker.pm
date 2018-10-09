=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::IdMapping::RNAProductPairMaker - produces pairs of
functionally similar mature RNA products from two ID-mapped transcripts

=head1 DESCRIPTION

Mature RNA products of a transcript are more complicated to perform
stable-ID mapping on than (canonical) translations because a
transcript can have multiple mature RNA products associated with
it. For example, in case of MicroRNAs there are typically two strands
originating from a single precursor hairpin - one from the 5' side and
one from the 3' side.

In order for ID mapping to work correctly on RNAProducts it is
necessary to classify those products by function before attempting to
match them between the source and the target transcript. This module
has been created in order to keep all classification logic in one
place instead of repeating it it in all parts of ID-mapping code that
require it.

=head1 SYNOPSIS

  use Bio::EnsEMBL::IdMapping::RNAProductPairMaker qw(pair_rnaproducts);

  my $rnaproduct_mappings = Bio::EnsEMBL::IdMapping::MappingList->new(...);
  ...
  my $i = pair_rnaproducts($transcript_mappings->get_all_Entries(), $object_cache, sub {
    my ($src_product, $tgt_product, $transcript_entry) = @_;

    my $rp_entry = Bio::EnsEMBL::IdMapping::Entry->new_fast([
      $src_product->id(), $tgt_product->id(), $transcript_entry->score()
    ]);
    $rnaproduct_mappings->add_Entry($rp_entry);

  });

=cut



package Bio::EnsEMBL::IdMapping::RNAProductPairMaker;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw( throw );


=head2 pair_rnaproducts

  Arg [1]    : ArrayRef $mapping_list List of transcript mappings to
               iterate over
  Arg [2]    : Bio::EnsEMBL::IdMapping::Cache $cache A cache object

  Arg [3]    : CodeRef $match_action The callback to be executed on
               a matching pair of rnaproducts. Will be passed three
               arguments: the source rnaproduct, the target rnaproduct,
               and the transcript mapping from which they have come.
  Example    : pair_rnaproducts($transcript_mappings->get_all_Entries(),
                                $self->cache,
                                sub {
                                  print "It's a match: " . join(' ', @_);
                                });
  Description: Iterate over the given list of transcript mappings and
               execute the callback for every valid source-target pair
               of rnaproducts. What constitutes valid depends on the
               specific type of mature RNA products in question, for
               instance for miRNAs a pair is valid if it comes from
               the same arm of the precursor hairpin in both the
               source and the target.
               Returns the number of transcript mappings with no
               rnaproducts on at least one (i.e. source or target)
               side. Note that this number might be smaller than the
               number of transcript mappings with no corresponding
               rnaproduct mappings.
  Returntype : Integer
  Exceptions : none
  Caller     : Bio::EnsEMBL::IdMapping::InternalIdMapper,
               Bio::EnsEMBL::IdMapping::StableIdMapper
  Status     : Stable

=cut

sub pair_rnaproducts {
  my ($mapping_list, $cache, $match_action) = @_;

  my $i = 0;
  foreach my $entry (@{ $mapping_list }) {

    my $source_rps = $cache->get_by_key('transcripts_by_id',
      'source', $entry->source)->get_all_RNAProducts();
    my $target_rps = $cache->get_by_key('transcripts_by_id',
      'target', $entry->target)->get_all_RNAProducts();

    if ((scalar @{$source_rps} != 0) && (scalar @{$target_rps} != 0)) {
      my %src_rp_map = _classify_rnaproducts('source', $source_rps);
      my %tgt_rp_map = _classify_rnaproducts('target', $target_rps);

      while (my ($rp_type, $type_submap) = each %src_rp_map) {

        if ($rp_type eq 'Bio::EnsEMBL::MicroRNA') {

          # Add a mapping for each existing pair of miRNA on the same arm
          # of the hairpin
          while (my ($mirna_arm, $src_product) = each %{ $type_submap }) {
            my $tgt_product = $src_rp_map{$rp_type}->{$mirna_arm};
            if ($src_product && $tgt_product) {
              &{$match_action}($src_product, $tgt_product, $entry);
            }
          }

        }

        # _classify_rnaproducts() will have filtered out unsupported types by now

      }
    }
    else {
      ++$i;
    }

  }

  return $i;
}


# _classify_rnaproducts

# Arg [1]    : String $message_tag A label to include in exception
#              messages to make them more informative
# Arg [2]    : ArrayRef $rps List of RNAProduct objects to process
# Description: PRIVATE For a given list of rnaproduct entries, tag
#              them with type and appropriate type-specific
#              properties. Prevents having to make multiple passes
#              through rnaproduct lists in the course of
#              (many-to-many) comparisons between rnaproducts of
#              source and target transcript, moreover it explicitly
#              rejects unsupported rnaproduct types and malformed
#              input.
# Returntype : Hash
# Exceptions : Throws on malformed data for known rnaproduct types,
#              and upon encountering an unsupported type.
# Caller     : internal
# Status     : Stable
sub _classify_rnaproducts {
  my ($message_tag, $rps) = @_;

  my %rp_map;
  foreach my $rnaproduct (@{ $rps }) {
    my $class_name = ref($rnaproduct);

    if ($class_name eq 'Bio::EnsEMBL::MicroRNA') {
      my $mirna_arm = $rnaproduct->arm();

      # Sanity checks
      if (!defined($mirna_arm)) {
        throw("Malformed data: ${message_tag} MicroRNA '"
              . $rnaproduct->display_id() . "' has no arm attribute");
      }
      if (exists($rp_map{$class_name}->{$mirna_arm})) {
        throw("Malformed data: ${message_tag} transcript has multiple "
              . "MicroRNA objects with arm=${mirna_arm}");
      }

      $rp_map{$class_name} = {
        $mirna_arm => $rnaproduct,
      };

    }
    else {
      throw("I do not know how to map RNA-product type ${class_name}");
    }
  }

  return %rp_map;
}


1;
