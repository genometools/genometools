/*
  Copyright (c) 2010-2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010-2013 Center for Bioinformatics, University of Hamburg

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

#ifndef LTRHARVEST_STREAM_H
#define LTRHARVEST_STREAM_H

#include "match/xdrop.h"
#include "extended/node_stream_api.h"
#include "ltr/ltr_four_char_motif.h"

/* implements the ``genome stream'' interface */
typedef struct GtLTRharvestStream GtLTRharvestStream;

const GtNodeStreamClass* gt_ltrharvest_stream_class(void);
GtNodeStream* gt_ltrharvest_stream_new(GtStr *str_indexname,
                                       GtRange searchrange,
                                       GtUword minseedlength,
                                       GtUword minltrlength,
                                       GtUword maxltrlength,
                                       GtUword mindistance,
                                       GtUword maxdistance,
                                       double similaritythreshold,
                                       int xdropbelowscore,
                                       GtXdropArbitraryscores arbitscores,
                                       GtLTRFourCharMotif *motif,
                                       bool verbosemode,
                                       bool nooverlaps,
                                       bool bestoverlaps,
                                       bool scanfile,
                                       GtUword offset,
                                       unsigned int minlengthTSD,
                                       unsigned int maxlengthTSD,
                                       GtUword vicinity,
                                       GtError *err);

const GtEncseq* gt_ltrharvest_stream_get_encseq(GtNodeStream *ltrh_stream);
void            gt_ltrharvest_stream_disable_seqids(GtLTRharvestStream*);
void            gt_ltrharvest_stream_enable_seqids(GtLTRharvestStream*);
void            gt_ltrharvest_stream_disable_md5_seqids(GtLTRharvestStream*);
void            gt_ltrharvest_stream_enable_md5_seqids(GtLTRharvestStream*);

#endif
