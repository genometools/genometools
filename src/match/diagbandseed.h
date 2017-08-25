/*
  Copyright (c) 2015-2016 Joerg Winkler <j.winkler@posteo.de>
  Copyright (c) 2016 Stefan Kurtz  <kurtz@zbh.uni-hamburg.de>
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
#include "match/seed_extend_parts.h"
#include "match/querymatch-display.h"
#include "match/xdrop.h"

typedef struct GtDiagbandseedInfo GtDiagbandseedInfo;
typedef struct GtDiagbandseedExtendParams GtDiagbandseedExtendParams;

typedef enum
{ /* keep the order consistent with gt_base_list_arguments */
  GT_DIAGBANDSEED_BASE_LIST_STRUCT,
  GT_DIAGBANDSEED_BASE_LIST_ULONG,
  GT_DIAGBANDSEED_BASE_LIST_BYTESTRING,
  GT_DIAGBANDSEED_BASE_LIST_UNDEFINED
} GtDiagbandseedBaseListType;

/* Run the whole algorithm. */
int gt_diagbandseed_run(const GtDiagbandseedInfo *arg,
                        const GtSequencePartsInfo *aseqranges,
                        const GtSequencePartsInfo *bseqranges,
                        const GtUwordPair *pick,
                        GtError *err);

/* The constructor for GtDiagbandseedInfo*/
GtDiagbandseedInfo *gt_diagbandseed_info_new(const GtEncseq *aencseq,
                                             const GtEncseq *bencseq,
                                             GtUword maxfreq,
                                             GtUword memlimit,
                                             unsigned int spacedseedweight,
                                             unsigned int seedlength,
                                             bool norev,
                                             bool nofwd,
                                             const GtRange *seedpairdistance,
                                             GtDiagbandseedBaseListType splt,
                                             GtDiagbandseedBaseListType kmplt,
                                             bool verify,
                                             bool verbose,
                                             bool debug_kmer,
                                             bool debug_seedpair,
                                             bool use_kmerfile,
                                             bool trimstat_on,
                                             GtUword maxmat,
                                             const GtStr *chainarguments,
                                             const GtStr
                                               *diagband_statistics_arg,
                                             size_t file_buffer_size,
                                             bool snd_pass,
                                             bool inseqseeds,
                                             const GtDiagbandseedExtendParams
                                               *extp);

const char *gt_diagbandseed_splt_comment(void);

const char *gt_diagbandseed_kmplt_comment(void);

GtDiagbandseedBaseListType gt_diagbandseed_base_list_get(
                                   bool with_splt,
                                   const char *base_list_string,GtError *err);

typedef struct GtAniAccumulate GtAniAccumulate;

GtAniAccumulate *gt_ani_accumulate_new(const GtEncseq *a_encseq,
                                       const GtEncseq *b_encseq,
                                       bool aseqfirstrun);

void gt_ani_accumulate_delete(GtAniAccumulate *ani_accumulate);

/* The constructor for GtDiagbandseedExtendParams*/
GtDiagbandseedExtendParams *gt_diagbandseed_extend_params_new(
                                GtUword userdefinedleastlength,
                                GtUword errorpercentage,
                                double evalue_threshold,
                                GtUword logdiagbandwidth,
                                GtUword mincoverage,
                                const GtSeedExtendDisplayFlag *out_display_flag,
                                GtUword use_apos,
                                GtXdropscore xdropbelowscore,
                                bool extendgreedy,
                                bool extendxdrop,
                                GtUword maxalignedlendifference,
                                GtUword history_size,
                                GtUword perc_mat_history,
                                GtExtendCharAccess a_extend_char_access,
                                GtExtendCharAccess b_extend_char_access,
                                bool cam_generic,
                                GtUword sensitivity,
                                double matchscore_bias,
                                bool weakends,
                                bool benchmark,
                                bool always_polished_ends,
                                bool verify_alignment,
                                bool only_selected_seqpairs,
                                GtAniAccumulate *ani_accumulate);

bool gt_diagbandseed_derive_maxmat_show(GtUword maxmat);

/* The destructors */
void gt_diagbandseed_info_delete(GtDiagbandseedInfo *info);

void gt_diagbandseed_extend_params_delete(GtDiagbandseedExtendParams *extp);

#endif
