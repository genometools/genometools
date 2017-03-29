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
#include "core/str_array_api.h"
#include "core/error_api.h"

typedef struct GtSeedExtendDisplayFlag GtSeedExtendDisplayFlag;

GtUword gt_querymatch_display_alignmentwidth(const GtSeedExtendDisplayFlag
                                                *display_flag);

void gt_querymatch_display_seedpos_relative_set(GtSeedExtendDisplayFlag
                                                *display_flag,
                                                bool a_is_rel,
                                                bool b_is_rel);

bool gt_querymatch_display_seedpos_a_relative(const GtSeedExtendDisplayFlag
                                                *display_flag);

bool gt_querymatch_display_seedpos_b_relative(const GtSeedExtendDisplayFlag
                                                *display_flag);

GtSeedExtendDisplayFlag *gt_querymatch_display_flag_new(const GtStrArray
                                                          *display_args,
                                                        GtError *err);

void gt_querymatch_display_flag_delete(GtSeedExtendDisplayFlag *display_flag);

const char *gt_querymatch_display_help(void);

void gt_querymatch_fields_approx_output(const GtSeedExtendDisplayFlag
                                         *display_flag,FILE *stream);

void gt_querymatch_fields_exact_output(FILE *stream);

bool gt_querymatch_alignment_display(const GtSeedExtendDisplayFlag
                                     *display_flag);

const unsigned int *gt_querymatch_display_order(unsigned int *numcolumns,
                                                const GtSeedExtendDisplayFlag
                                                   *display_flag);

#include "match/se-display-fwd.inc"

#endif
