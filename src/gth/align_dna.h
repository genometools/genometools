/*
  Copyright (c) 2003-2010 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef ALIGN_DNA_H
#define ALIGN_DNA_H

#include "gth/align_common.h"
#include "gth/dp_options_core.h"
#include "gth/dp_options_est.h"
#include "gth/dp_options_postpro.h"
#include "gth/path_matrix.h"
#include "gth/sa.h"

#define DASH_DASH_WEIGHT        0.0

/* XXX: precompute this for performance reasons (part of core DP) */
#define ADDOUTPUTWEIGHT(VAR,N,M)\
        if ((N) < (gen_alphabet_mapsize-1))\
        {\
          if ((M) < (gen_alphabet_mapsize-1))\
          {\
            if ((N)==(M))\
            {\
              VAR += dp_options_est->identityweight;\
            }\
            else\
            {\
              VAR += dp_options_est->mismatchweight;\
            }\
          }\
          else if ((M) == DASH)\
          {\
            VAR += dp_options_est->deletionweight;\
          }\
          else\
          {\
            VAR += dp_options_est->undetcharweight;\
          }\
        }\
        else if ((N) == DASH)\
        {\
          if ((M) == DASH)\
          {\
            /* N=M=DASH: can happen in gthcomputescores() */\
            VAR += DASH_DASH_WEIGHT;\
          }\
          else\
          {\
            VAR += dp_options_est->deletionweight;\
          }\
        }\
        else\
        {\
          if ((M) == DASH)\
          {\
            VAR += dp_options_est->deletionweight;\
          }\
          else\
          {\
            VAR += dp_options_est->undetcharweight;\
          }\
        }

typedef struct GthDPMatrix GthDPMatrix;

typedef void (*GthDNACompletePathMatrixJT)(GthDPMatrix *dpm,
                                           const unsigned char *gen_seq_tran,
                                           const unsigned char *ref_seq_tran,
                                           GtUword genomic_offset,
                                           GtAlphabet *gen_alphabet,
                                           GthDPParam *dp_param,
                                           GthDPOptionsEST *dp_options_est,
                                           GthDPOptionsCore *dp_options_core,
                                           GthJumpTable *jumpt_table,
                                           GtArray *gen_ranges,
                                           GtUword ref_dp_length,
                                           GtUword ref_offset,
                                           GthPathMatrix **pm);
#define ADDOUTPUTWEIGHTIDENTITY(VAR,N)\
        if ((N) < (gen_alphabet_mapsize-1))\
        {\
          VAR += dp_options_est->identityweight;\
        }\
        else if ((N) == DASH)\
        {\
          VAR += DASH_DASH_WEIGHT;\
        }\
        else\
        {\
          VAR += dp_options_est->undetcharweight;\
        }

/* The following function implements the Spliced Alignment of Genomic DNA with
   cDNA, as described by Usuka, Zhu and Brendel. */
int gth_align_dna(GthSA*,
                  GtArray *gen_ranges,
                  const unsigned char *gen_seq_tran,
                  const unsigned char *gen_seq_orig,
                  const unsigned char *ref_seq_tran,
                  const unsigned char *ref_seq_orig,
                  GtUword ref_dp_length,
                  GtAlphabet *gen_alphabet,
                  GtAlphabet *ref_alphabet,
                  bool introncutout,
                  GtUword autoicmaxmatrixsize,
                  bool showeops,
                  bool comments,
                  bool gs2out,
                  const GtRange *gen_seq_bounds,
                  GthSpliceSiteModel *splice_site_model,
                  GthDPOptionsCore *dp_options_core,
                  GthDPOptionsEST *dp_options_est,
                  GthDPOptionsPostpro *dp_options_postpro,
                  GthDNACompletePathMatrixJT complete_path_matrix_jt,
                  GthJumpTable *jump_table,
                  GtUword ref_offset,
                  GthStat*,
                  GtFile*);

/* can return NULL */
GthSA* gth_align_dna_simple(GthInput *input,
                            const GtRange *gen_range,
                            GtUword gen_file_num,
                            GtUword gen_seq_num,
                            bool gen_strand_forward,
                            GtUword ref_file_num,
                            GtUword ref_seq_num,
                            GthSpliceSiteModel *splice_site_model);

void gth_show_backtrace_matrix(GthPath **path,
                               GtUword gen_dp_length,
                               GtUword ref_dp_length,
                               const GtRange *btmatrixgenrange,
                               const GtRange *btmatrixrefrange);

#endif
