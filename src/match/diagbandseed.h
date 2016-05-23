/*
  Copyright (c) 2015-2016 Joerg Winkler <j.winkler@posteo.de>
  Copyright (c) 2015-2016 Center for Bioinformatics, University of Hamburg

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

#ifndef DIAGBANDSEED_H
#define DIAGBANDSEED_H
#include <stdbool.h>
#include "core/encseq_api.h"
#include "core/error_api.h"
#include "core/range_api.h"
#include "core/types_api.h"
#include "match/ft-front-prune.h"
#include "match/xdrop.h"

typedef struct GtDiagbandseedInfo GtDiagbandseedInfo;
typedef struct GtDiagbandseedExtendParams GtDiagbandseedExtendParams;

/* Run the whole algorithm. */
int gt_diagbandseed_run(const GtDiagbandseedInfo *arg,
                        const GtRange *aseqranges,
                        const GtRange *bseqranges,
                        GtError *err);

/* The constructor for GtDiagbandseedInfo*/
GtDiagbandseedInfo *gt_diagbandseed_info_new(GtEncseq *aencseq,
                                             GtEncseq *bencseq,
                                             GtUword maxfreq,
                                             GtUword memlimit,
                                             unsigned int seedlength,
                                             bool norev,
                                             bool nofwd,
                                             GtRange *seedpairdistance,
                                             bool verify,
                                             bool verbose,
                                             bool debug_kmer,
                                             bool debug_seedpair,
                                             bool use_kmerfile,
                                             GtDiagbandseedExtendParams *extp,
                                             GtUword anumseqranges,
                                             GtUword bnumseqranges);

/* The constructor for GtDiagbandseedExtendParams*/
GtDiagbandseedExtendParams *gt_diagbandseed_extend_params_new(
                              GtUword errorpercentage,
                              GtUword userdefinedleastlength,
                              GtUword logdiagbandwidth,
                              GtUword mincoverage,
                              unsigned int display_flag,
                              bool use_apos,
                              GtXdropscore xdropbelowscore,
                              bool extendgreedy,
                              bool extendxdrop,
                              GtUword maxalignedlendifference,
                              GtUword history_size,
                              GtUword perc_mat_history,
                              GtExtendCharAccess extend_char_access,
                              GtUword sensitivity,
                              double matchscore_bias,
                              bool weakends,
                              bool benchmark,
                              GtUword alignmentwidth,
                              bool always_polished_ends);

/* The destructors */
void gt_diagbandseed_info_delete(GtDiagbandseedInfo *info);

void gt_diagbandseed_extend_params_delete(GtDiagbandseedExtendParams *extp);
#endif
