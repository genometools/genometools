/*
  Copyright (c) 2008-2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2009 Center for Bioinformatics, University of Hamburg

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

#ifndef LTRDIGEST_STREAM_H
#define LTRDIGEST_STREAM_H

#include "core/error_api.h"
#include "extended/node_stream_api.h"
#include "ltr/pbs.h"
#include "ltr/ppt.h"
#ifdef HAVE_HMMER
#include "ltr/pdom.h"
#endif
#include "core/encseq.h"

/* implements the ``node stream'' interface */
typedef struct GtLTRdigestStream GtLTRdigestStream;

const GtNodeStreamClass* gt_ltrdigest_stream_class(void);

GtNodeStream* gt_ltrdigest_stream_new(GtNodeStream *in_stream,
                                      int tests_to_run,
                                      GtEncseq *encseq,
                                      GtPBSOptions *pbs_opts,
                                      GtPPTOptions *ppt_opts,
#ifdef HAVE_HMMER
                                      GtPdomOptions *pdom_opts,
#endif
                                      GtError *err);

#endif
