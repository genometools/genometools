/*
  Copyright (c) 2003-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef COMPUTE_SCORES_H
#define COMPUTE_SCORES_H

#include <stdbool.h>
#include "core/trans_table_api.h"
#include "gth/align_common.h"
#include "gth/dp_scores_protein.h"
#include "gth/sa.h"
#include "gth/spliced_seq.h"

void gth_compute_scores(GthSA *sa,
                        bool proteineop,
                        GthDPParam *dp_param,
                        void *dp_options_est,
                        const unsigned char *gen_seq_tran,
                        const unsigned char *ref_seq_tran,
                        const unsigned char *ref_seq_orig,
                        const GtTransTable *transtable,
                        unsigned long gen_dp_start,
                        unsigned long scoreminexonlen,
                        bool introncutout,
                        bool gs2out,
                        GthSplicedSeq *spliced_seq,
                        unsigned long ref_dp_length,
                        GtAlphabet *gen_alphabet,
                        GtAlphabet *ref_alphabet,
                        GthDPScoresProtein *dp_scores_protein);

#endif
