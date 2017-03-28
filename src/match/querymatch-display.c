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

typedef struct
{
  const char *name;
  int rank;
  bool incolumn;
} GtSEdisplayStruct;

struct GtSeedExtendDisplayFlag
{
  unsigned int flags,
               nextfree,
               order[GT_DISPLAY_LARGEST_FLAG+1];
  bool a_seedpos_relative, b_seedpos_relative;
  GtUword alignmentwidth;
};

static unsigned int gt_display_mask(int shift)
{
  gt_assert(shift <= GT_DISPLAY_LARGEST_FLAG);
  return 1U << shift;
}

static bool gt_querymatch_display_on(const GtSeedExtendDisplayFlag
                                       *display_flag,
                                     GtSeedExtendDisplay_enum display)
{
  return (display_flag != NULL &&
          (display_flag->flags & gt_display_mask(display))) ? true : false;
}

#include "match/se-display.inc"

static const GtSEdisplayStruct *gt_display_arg_get(const char *str,
                                                   size_t cmplen)
{
  const GtSEdisplayStruct *left = gt_display_arguments_table,
                          *right = gt_display_arguments_table +
                                   sizeof gt_display_arguments_table/
                                   sizeof gt_display_arguments_table[0] - 1;

  while (left <= right)
  {
    const GtSEdisplayStruct *mid = left + (right - left + 1)/2;
    int cmp = cmplen == 0 ? strcmp(str,mid->name)
                          : strncmp(str,mid->name,cmplen);

    if (cmp < 0)
    {
      right = mid - 1;
    } else
    {
      if (cmp > 0)
      {
        left = mid + 1;
      } else
      {
        return mid;
      }
    }
  }
  return NULL;
}

static void gt_querymatch_display_flag_add(GtSeedExtendDisplayFlag
                                            *display_flag,unsigned int flag)
{
  display_flag->flags |= gt_display_mask(flag);
  display_flag->order[display_flag->nextfree++] = flag;
}

GtSeedExtendDisplayFlag *gt_querymatch_display_flag_new(void)
{
  GtSeedExtendDisplayFlag *display_flag = gt_malloc(sizeof *display_flag);

  display_flag->alignmentwidth = 0;
  display_flag->a_seedpos_relative = true; /* as bytes is default access mode */
  display_flag->b_seedpos_relative = true;
  display_flag->nextfree = 0;
  display_flag->flags = 0;
  gt_querymatch_display_flag_add(display_flag,Gt_S_len_display);
  gt_querymatch_display_flag_add(display_flag,Gt_S_seqnum_display);
  gt_querymatch_display_flag_add(display_flag,Gt_S_start_display);
  gt_querymatch_display_flag_add(display_flag,Gt_Strand_display);
  gt_querymatch_display_flag_add(display_flag,Gt_Q_len_display);
  gt_querymatch_display_flag_add(display_flag,Gt_Q_seqnum_display);
  gt_querymatch_display_flag_add(display_flag,Gt_Q_start_display);
  gt_querymatch_display_flag_add(display_flag,Gt_Score_display);
  gt_querymatch_display_flag_add(display_flag,Gt_Editdist_display);
  gt_querymatch_display_flag_add(display_flag,Gt_Identity_display);
  return display_flag;
}

GtUword gt_querymatch_display_alignmentwidth(const GtSeedExtendDisplayFlag
                                                *display_flag)
{
  return (display_flag == NULL) ? 0 : display_flag->alignmentwidth;
}

bool gt_querymatch_alignment_display(const GtSeedExtendDisplayFlag
                                     *display_flag)
{
  return (display_flag != NULL &&
          display_flag->alignmentwidth > 0) ? true : false;
}

void gt_querymatch_display_seedpos_relative_set(GtSeedExtendDisplayFlag
                                                *display_flag,
                                                bool a_is_rel,
                                                bool b_is_rel)
{
  gt_assert(display_flag != NULL);
  display_flag->a_seedpos_relative = a_is_rel;
  display_flag->b_seedpos_relative = b_is_rel;
}

bool gt_querymatch_display_seedpos_a_relative(const GtSeedExtendDisplayFlag
                                                *display_flag)
{
  return (display_flag != NULL && display_flag->a_seedpos_relative) ? true
                                                                    : false;
}

bool gt_querymatch_display_seedpos_b_relative(const GtSeedExtendDisplayFlag
                                                *display_flag)
{
  return (display_flag != NULL && display_flag->b_seedpos_relative) ? true
                                                                    : false;
}

void gt_querymatch_fields_approx_output(const GtSeedExtendDisplayFlag
                                         *display_flag,FILE *stream)
{
  GtStr *add_column_header;

  fprintf(stream,"# Fields: ");
  add_column_header = gt_querymatch_column_header(display_flag);
  gt_assert(gt_str_length(add_column_header) > 0);
  fputs(gt_str_get(add_column_header),stream);
  fputc('\n',stream);
  gt_str_delete(add_column_header);
}

void gt_querymatch_fields_exact_output(FILE *stream)
{
  fprintf(stream,"# Fields: s.len, s.seqnum, s.start, strand, q.seqnum, "
                 "q.start\n");
}

static int gt_querymatch_display_flag_set(GtWord *parameter,
                                          GtSeedExtendDisplayFlag *display_flag,
                                          const char *arg,
                                          GtError *err)
{
  const char *exclude_list[] = {"alignment","cigar"};
  size_t ex_idx, numexcl = sizeof exclude_list/sizeof exclude_list[0];
  const GtSEdisplayStruct *dstruct;
  const char *ptr;
  size_t cmplen = 0;

  gt_assert(display_flag != NULL && numexcl % 2 == 0);
  ptr = strchr(arg,'=');
  if (ptr != NULL)
  {
    cmplen = (size_t) (ptr - arg);
    if (sscanf(ptr+1,GT_WD,parameter) != 1)
    {
      gt_error_set(err,"illegal argument \"%s\" to option -outfmt: "
                       "expect integer following symbol =",arg);
      return -1;
    }
  }
  dstruct = gt_display_arg_get(arg,cmplen);
  if (dstruct == NULL)
  {
    gt_error_set(err,"illegal identifier \"%s\" as argument of options "
                     "-outfmt, possible identifiers are: %s",arg,
                     GT_SE_POSSIBLE_DISPLAY_ARGS);
    return -1;
  }
  display_flag->flags |= gt_display_mask((int) dstruct->rank);
  for (ex_idx = 0; ex_idx < numexcl; ex_idx+=2)
  {
    const GtSEdisplayStruct
      *dstruct0 = gt_display_arg_get(exclude_list[ex_idx],0),
      *dstruct1 = gt_display_arg_get(exclude_list[ex_idx+1],0);

    gt_assert(dstruct0 != NULL && dstruct1 != NULL);
    if ((display_flag->flags & gt_display_mask(dstruct0->rank)) &&
        (display_flag->flags & gt_display_mask(dstruct1->rank)))
    {
      gt_error_set(err,"argument \"%s\" and \"%s\" of option -outfmt exclude "
                       "each other",exclude_list[ex_idx],
                                    exclude_list[ex_idx+1]);
      return -1;
    }
  }
  return cmplen > 0 ? 1 : 0;
}

int gt_querymatch_display_flag_args_set(GtSeedExtendDisplayFlag *display_flag,
                                        const GtStrArray *display_args,
                                        GtError *err)
{
  bool haserr = false;
  GtUword da_idx;

  gt_assert(display_flag != NULL);
  for (da_idx = 0; da_idx < gt_str_array_size(display_args); da_idx++)
  {
    const char *da = gt_str_array_get(display_args,da_idx);
    GtWord parameter;
    int ret = gt_querymatch_display_flag_set(&parameter,display_flag,da,err);
    switch (ret)
    {
      case 0:
        break;
      case 1:
        /* the only flag with a parameter is Gt_Alignment_display */
        if (parameter < 0)
        {
          gt_error_set(err,"integer following \"alignment=\" must be positive");
          haserr = true;
        } else
        {
          gt_assert(display_flag->flags &
                    gt_display_mask(Gt_Alignment_display));
          display_flag->alignmentwidth = (GtUword) parameter;
        }
        break;
      default:
        haserr = true;
    }
  }
  if (!haserr &&
      (display_flag->flags & gt_display_mask(Gt_Alignment_display)) &&
      display_flag->alignmentwidth == 0)
  {
    display_flag->alignmentwidth = 60; /* this is the default alignment width */
  }
  return haserr ? -1 : 0;
}

void gt_querymatch_display_flag_delete(GtSeedExtendDisplayFlag *display_flag)
{
  if (display_flag != NULL)
  {
    gt_free(display_flag);
  }
}
