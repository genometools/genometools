/*
  Copyright (c) 2017 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2017 Center for Bioinformatics, University of Hamburg

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
#include <stdbool.h>
#include <string.h>
#include "core/ma_api.h"
#include "core/assert_api.h"
#include "match/querymatch-display.h"

typedef enum
{
  Gt_Seed_display,
  Gt_Seqlength_display,
  Gt_Evalue_display,
  Gt_Seqdesc_display,
  Gt_Bitscore_display
} GtSeedExtendDisplay_enum;

struct GtSeedExtendDisplayFlag
{
  unsigned int flags;
  GtUword alignmentwidth;
};

static bool gt_querymatch_display_on(const GtSeedExtendDisplayFlag
                                       *display_flag,
                                     GtSeedExtendDisplay_enum display)
{
  gt_assert((int) display <= Gt_Bitscore_display);
  return (display_flag != NULL &&
          (display_flag->flags & (1U << (int) display))) ? true : false;
}

bool gt_querymatch_seed_display(const GtSeedExtendDisplayFlag *display_flag)
{
  return gt_querymatch_display_on(display_flag,Gt_Seed_display);
}

bool gt_querymatch_evalue_display(const GtSeedExtendDisplayFlag *display_flag)
{
  return gt_querymatch_display_on(display_flag,Gt_Evalue_display);
}

bool gt_querymatch_bit_score_display(const GtSeedExtendDisplayFlag
                                       *display_flag)
{
  return gt_querymatch_display_on(display_flag,Gt_Bitscore_display);
}

bool gt_querymatch_seq_desc_display(const GtSeedExtendDisplayFlag *display_flag)
{
  return gt_querymatch_display_on(display_flag,Gt_Seqdesc_display);
}

bool gt_querymatch_seqlength_display(const GtSeedExtendDisplayFlag
                                       *display_flag)
{
  return gt_querymatch_display_on(display_flag,Gt_Seqlength_display);
}

GtUword gt_querymatch_display_alignmentwidth(const GtSeedExtendDisplayFlag
                                                *display_flag)
{
  return (display_flag == NULL) ? 0 : display_flag->alignmentwidth;
}

GtStr *gt_querymatch_column_header(const GtSeedExtendDisplayFlag *display_flag)
{
  GtStr *str = gt_str_new();

  if (gt_querymatch_seqlength_display(display_flag))
  {
    gt_str_append_cstr(str," aseqlen bseqlen");
  }
  if (gt_querymatch_evalue_display(display_flag))
  {
    gt_str_append_cstr(str," evalue");
  }
  if (gt_querymatch_bit_score_display(display_flag))
  {
    gt_str_append_cstr(str," bit-score");
  }
  return str;
}

const char *gt_querymatch_display_help(void)
{
  return "specify what additional values in matches are displayed\n"
         "seed:      display the seed of the match\n"
         "seqlength: display length of sequences in which\n"
         "           the two match-instances occur\n"
         "evalue:    display evalue\n"
         "seq-desc:  display sequence description instead of numbers\n"
         "bit-score: display bit score";
}

static bool gt_querymatch_display_flag_set(GtSeedExtendDisplayFlag
                                             *display_flag,
                                           const char *arg)
{
  const char *display_strings[]
    = {"seed","seqlength","evalue","seq-desc","bit-score"};
  size_t ds_idx, numofds = sizeof display_strings/sizeof display_strings[0];
  bool found = false;

  gt_assert(display_flag != NULL &&
            numofds == (size_t) Gt_Bitscore_display + 1);
  for (ds_idx = 0; ds_idx < numofds; ds_idx++)
  {
    if (strcmp(arg,display_strings[ds_idx]) == 0)
    {
      display_flag->flags |= (1U << ds_idx);
      found = true;
    }
  }
  return found;
}

void gt_querymatch_display_alignmentwidth_set(GtSeedExtendDisplayFlag
                                                  *display_flag,
                                              GtUword alignmentwidth)
{
  gt_assert(display_flag != NULL);
  display_flag->alignmentwidth = alignmentwidth;
}

GtSeedExtendDisplayFlag *gt_querymatch_display_flag_new(void)
{
  GtSeedExtendDisplayFlag *display_flag = gt_malloc(sizeof *display_flag);

  display_flag->flags = 0;
  display_flag->alignmentwidth = 0;
  return display_flag;
}

int gt_querymatch_display_flag_args_set(
               GtSeedExtendDisplayFlag *display_flag,
               const GtStrArray *display_args,
               GtError *err)
{
  GtUword da_idx;

  gt_assert(display_flag != NULL);
  for (da_idx = 0; da_idx < gt_str_array_size(display_args); da_idx++)
  {
    const char *da = gt_str_array_get(display_args,da_idx);

    if (!gt_querymatch_display_flag_set(display_flag,da))
    {
      gt_error_set(err,"illegal argument %s to option -display: "
                       " possible values are "
                       "seed, seqlength, evalue, seq-desc or bit-score",da);
      gt_free(display_flag);
      return -1;
    }
  }
  return 0;
}

void gt_querymatch_display_flag_delete(GtSeedExtendDisplayFlag *display_flag)
{
  if (display_flag != NULL)
  {
    gt_free(display_flag);
  }
}
