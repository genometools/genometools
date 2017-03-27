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

#include "match/se-display.inc"

const GtSEdisplayStruct *gt_display_arg_get(const char *str)
{
  const GtSEdisplayStruct *left = gt_display_arguments_table,
                          *right = gt_display_arguments_table +
                                   sizeof gt_display_arguments_table/
                                   sizeof gt_display_arguments_table[0] - 1;

  while (left <= right)
  {
    const GtSEdisplayStruct *mid = left + (right - left + 1)/2;
    int cmp = strcmp(str,mid->name);

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

GtSeedExtendDisplayFlag *gt_querymatch_display_flag_new(void)
{
  GtSeedExtendDisplayFlag *display_flag = gt_malloc(sizeof *display_flag);
  GtStrArray *default_display_args = gt_str_array_new();
  int ret;

  display_flag->flags = 0;
  display_flag->alignmentwidth = 0;
  display_flag->a_seedpos_relative = true; /* as bytes is default access mode */
  display_flag->b_seedpos_relative = true;
  gt_str_array_add_cstr(default_display_args,"s.len");
  gt_str_array_add_cstr(default_display_args,"s.seqnum");
  gt_str_array_add_cstr(default_display_args,"s.start");
  gt_str_array_add_cstr(default_display_args,"strand");
  gt_str_array_add_cstr(default_display_args,"q.len");
  gt_str_array_add_cstr(default_display_args,"q.seqnum");
  gt_str_array_add_cstr(default_display_args,"q.start");
  gt_str_array_add_cstr(default_display_args,"score");
  gt_str_array_add_cstr(default_display_args,"editdist");
  gt_str_array_add_cstr(default_display_args,"identity");
  ret = gt_querymatch_display_flag_args_set(display_flag,default_display_args,
                                            NULL);
  gt_assert(ret == 0);
  gt_str_array_delete(default_display_args);
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

static GtStrArray *gt_se_help2display_strings(const char *helpline)
{
  const char *ptr = helpline, *laststart = NULL;
  char findchar = '\n';
  GtStrArray *display_args = gt_str_array_new();

  while (true)
  {
    ptr = strchr(ptr,findchar);
    gt_assert(ptr != NULL);
    if (findchar == '\n')
    {
      if (*(ptr+1) == '\0')
      {
        break;
      }
      ptr++;
      if (*ptr != ' ')
      {
        laststart = ptr;
        findchar = ':';
      }
    } else
    {
      GtUword len = (GtUword) (ptr - laststart);
      gt_str_array_add_cstr_nt(display_args,laststart,len);
      findchar = '\n';
    }
  }
  return display_args;
}

static int gt_querymatch_display_flag_set(GtWord *parameter,
                                          const GtStrArray *display_strings,
                                          GtSeedExtendDisplayFlag *display_flag,
                                          const char *arg,
                                          GtError *err)
{
  GtUword ds_idx, numofds = gt_str_array_size(display_strings);
  const GtSeedExtendDisplay_enum exclude_list[] = {Gt_Alignment_display,
                                                   Gt_Cigar_display};
  size_t ex_idx, numexcl = sizeof exclude_list/sizeof exclude_list[0];
  bool identifier_okay = false, parameter_found = false;
  const char *ptr;
  size_t cmplen;

  gt_assert(display_flag != NULL &&
            numofds == (size_t) Gt_Bitscore_display + 1 &&
            numexcl % 2 == 0);
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
    parameter_found = true;
  } else
  {
    cmplen = 0;
  }
  for (ds_idx = 0; ds_idx < numofds; ds_idx++)
  {
    const char *dstring = gt_str_array_get(display_strings,ds_idx);
    int ret = (cmplen > 0) ? strncmp(arg,dstring,cmplen)
                           : strcmp(arg,dstring);
    if (ret == 0)
    {
      display_flag->flags |= gt_display_mask((int) ds_idx);
      identifier_okay = true;
      break;
    }
  }
  if (!identifier_okay)
  {
    GtStr *possible_args = gt_str_new();

    gt_str_append_cstr(possible_args,
                       "illegal identifier %s as argument of option -outfmt: "
                       "possible idenfifiers are: ");
    for (ds_idx = 0; ds_idx < numofds; ds_idx++)
    {
      const char *dstring = gt_str_array_get(display_strings,ds_idx);
      gt_str_append_cstr(possible_args,dstring);
      if (ds_idx < numofds - 1)
      {
        gt_str_append_cstr(possible_args,", ");
      }
    }
    gt_error_set(err,"illegal identifier %s as argument of options -outfmt,"
                     "possible identifiers are %s",arg,
                     gt_str_get(possible_args));
    gt_str_delete(possible_args);
    return -1;
  }
  for (ex_idx = 0; ex_idx < numexcl; ex_idx+=2)
  {
    if ((display_flag->flags & gt_display_mask(exclude_list[ex_idx])) &&
        (display_flag->flags & gt_display_mask(exclude_list[ex_idx+1])))
    {
      const char *d1 = gt_str_array_get(display_strings,exclude_list[ex_idx]);
      const char *d2 = gt_str_array_get(display_strings,exclude_list[ex_idx+1]);
      gt_error_set(err,"argument \"%s\" and \"%s\" of option -outfmt exclude "
                       "each other",d1,d2);
      return -1;
    }
  }
  return parameter_found ? 1 : 0;
}

int gt_querymatch_display_flag_args_set(GtSeedExtendDisplayFlag *display_flag,
                                        const GtStrArray *display_args,
                                        GtError *err)
{
  bool haserr = false;
  GtUword da_idx;
  GtStrArray *display_strings
    = gt_se_help2display_strings(gt_querymatch_display_help());

  gt_assert(display_flag != NULL);
  for (da_idx = 0; da_idx < gt_str_array_size(display_args); da_idx++)
  {
    const char *da = gt_str_array_get(display_args,da_idx);
    GtWord parameter;
    int ret = gt_querymatch_display_flag_set(&parameter,display_strings,
                                             display_flag,da,err);
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
  gt_str_array_delete(display_strings);
  return haserr ? -1 : 0;
}

void gt_querymatch_display_flag_delete(GtSeedExtendDisplayFlag *display_flag)
{
  if (display_flag != NULL)
  {
    gt_free(display_flag);
  }
}
