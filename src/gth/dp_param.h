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

#ifndef DP_PARAM_H
#define DP_PARAM_H

#include "core/alphabet.h"
#include "gth/gthchain.h"
#include "gth/splice_site_model.h"

typedef struct {
  GthFlt *log_Pdonor,                 /* donor site */
         *log_1minusPdonor,
         *log_Pacceptor,              /* acceptor site */
         *log_1minusPacceptor;
} GthDPParam;

/* Can return NULL */
GthDPParam* gth_dp_param_new(GtArray *ranges,
                             const unsigned char *gen_seq_tran,
                             const GtRange *gen_seq_bounds,
                             GthSpliceSiteModel*,
                             GtAlphabet *gen_alphabet);
/* Can return NULL */
GthDPParam* gth_dp_param_new_with_range(unsigned long left,
                                        unsigned long right,
                                        const unsigned char *gen_seq_tran,
                                        const GtRange *gen_seq_bounds,
                                        GthSpliceSiteModel*,
                                        GtAlphabet *gen_alphabet);
void gth_dp_param_delete(GthDPParam *dp_param);

#endif
