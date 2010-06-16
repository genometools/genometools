/*
  Copyright (c) 2004-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2004-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef ALIGN_PROTEIN_H
#define ALIGN_PROTEIN_H

#include "core/trans_table_api.h"
#include "gth/dp_options_core.h"
#include "gth/dp_options_postpro.h"
#include "gth/spliced_seq.h"
#include "gth/align_common.h"
#include "gth/dp_scores_protein.h"
#include "gth/sa.h"

#define PROTEIN_NUMOFSCORETABLES  4

/* the following type bundles nearly all sahmtp input variables, except
   for gen_seq_tran and gen_dp_length */
typedef struct {
  const unsigned char *ref_seq_orig; /* pointer to original
                                                     reference sequence */
  GtScoreMatrix *score_matrix;      /* the amino acid substitution matrix */
  GtAlphabet *score_matrix_alpha;   /* alphabet used for the scoring matrix */
} GthAlignInputProtein;

unsigned char gthgetcodon(unsigned char genomicchar1,
                          unsigned char genomicchar2,
                          unsigned char genomicchar3,
                          const GtUchar *gen_alphabet_characters,
                          const GtTransTable *transtable);

/* The following function implements the Spliced Alignment of Genomic DNA with
   protein, as described by Usuka and Brendel. */
int gth_align_protein(GthSA*,
                      GtArray *gen_ranges,
                      const unsigned char *gen_seq_tran,
                      const unsigned char *ref_seq_tran,
                      const unsigned char *ref_seq_orig,
                      unsigned long referencelength,
                      GtAlphabet *gen_alphabet,
                      GtAlphabet *ref_alphabet,
                      GthInput *gth_input,
                      bool introncutout,
                      unsigned long autoicmaxmatrixsize,
                      bool proteinexonpenal,
                      bool showeops,
                      bool comments,
                      bool gs2out,
                      unsigned long translationtable,
                      const GtRange *gen_seq_bounds,
                      GthSpliceSiteModel *splice_site_model,
                      GthDPOptionsCore *dp_options_core,
                      GthDPOptionsPostpro *dp_options_postpro,
                      GthStat*,
                      GtFile*);

#endif
