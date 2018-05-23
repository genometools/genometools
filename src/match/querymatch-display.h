/*
  Copyright (c) 2007-2017 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2017 Center for Bioinformatics, University of Hamburg

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

#ifndef QUERYMATCH_DISPLAY_H
#define QUERYMATCH_DISPLAY_H
#include <stdbool.h>
#include "core/str_array_api.h"
#include "core/error_api.h"

typedef struct GtSeedExtendDisplayFlag GtSeedExtendDisplayFlag;

typedef enum
{
  GT_SEED_EXTEND_DISPLAY_SET_NO,
  GT_SEED_EXTEND_DISPLAY_SET_STANDARD,
  GT_SEED_EXTEND_DISPLAY_SET_EXACT
} GtSeedExtendDisplaySetMode;

GtSeedExtendDisplayFlag *gt_querymatch_display_flag_new(
                           const GtStrArray *display_args,
                           GtSeedExtendDisplaySetMode setmode,
                           GtError *err);

void gt_querymatch_display_flag_delete(GtSeedExtendDisplayFlag *display_flag);

void gt_querymatch_Fields_output(FILE *stream,
                                 const GtSeedExtendDisplayFlag *display_flag);

void gt_querymatch_Options_output(FILE *stream,int argc,const char **argv,
                                  bool idhistout,GtUword minidentity,
                                  GtUword historysize);

const unsigned int *gt_querymatch_display_order(GtUword *numcolumns,
                                                const GtSeedExtendDisplayFlag
                                                   *display_flag);

const char *gt_querymatch_display_help(void);

bool gt_querymatch_alignment_display(const GtSeedExtendDisplayFlag *);

bool gt_querymatch_run_aligner(const GtSeedExtendDisplayFlag *display_flag);

GtUword gt_querymatch_display_alignmentwidth(const GtSeedExtendDisplayFlag *);

GtUword gt_querymatch_trace_delta_display(const GtSeedExtendDisplayFlag *);

#include "match/se-display-fwd.inc"

const char *gt_querymatch_flag2name(GtSeedExtendDisplay_enum flag);

bool gt_querymatch_has_seed(const GtSeedExtendDisplayFlag *display_flag);

GtStrArray *gt_querymatch_read_Fields_line(const char *line_ptr);

#endif
