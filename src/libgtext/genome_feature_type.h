/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef GENOME_FEATURE_TYPE_H
#define GENOME_FEATURE_TYPE_H

/*
  Keep in sync with `genome_feature_type_strings' in `genome_feature_type.c'!
  The feature types have to be sorted. We assume ASCII encoding.
*/

typedef enum {
  gft_CDS,
  gft_EST_match,
  gft_LTR_retrotransposon,
  gft_SNP,
  gft_TF_binding_site,
  gft_cDNA_match,
  gft_exon,
  gft_five_prime_UTR,
  gft_five_prime_splice_site,
  gft_gene,
  gft_intron,
  gft_inverted_repeat,
  gft_long_terminal_repeat,
  gft_mRNA,
  gft_protein_match,
  gft_repeat_region,
  gft_target_site_duplication,
  gft_three_prime_UTR,
  gft_three_prime_splice_site,
  gft_transcript,
  undefined /* to have an ``undef'' enum, not convertable to and from cstr! */
} GenomeFeatureType;

/* Determine a genome feature type <gft> from the string <gft_string>.
   If such a feature does not exits, -1 is returned. */
int           genome_feature_type_get(GenomeFeatureType *gft,
                                      const char *gft_string);
const char*   genome_feature_type_get_cstr(GenomeFeatureType);
unsigned long genome_feature_type_num_of_features(void);

#endif
