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
#include <ctype.h>
#include <inttypes.h>
#include "core/ma_api.h"
#include "core/assert_api.h"
#include "match/querymatch-display.h"

typedef struct
{
  const char *name;
  GtSeedExtendDisplay_enum flag;
  bool incolumn;
} GtSEdisplayStruct;

struct GtSeedExtendDisplayFlag
{
  uint64_t flags;
  unsigned int order[GT_DISPLAY_LARGEST_FLAG+1];
  GtUword alignmentwidth, trace_delta, nextfree;
};

static uint64_t gt_display_mask(GtSeedExtendDisplay_enum flag)
{
  gt_assert(flag <= GT_DISPLAY_LARGEST_FLAG);
  return ((uint64_t) 1)  << flag;
}

static bool gt_querymatch_display_on(const GtSeedExtendDisplayFlag
                                       *display_flag,
                                     GtSeedExtendDisplay_enum display)
{
  return (display_flag != NULL &&
          (display_flag->flags & gt_display_mask(display))) ? true : false;
}

#include "match/se-display.inc"

bool gt_querymatch_has_seed(const GtSeedExtendDisplayFlag *display_flag)
{
  return gt_querymatch_seed_s_display(display_flag) &&
         gt_querymatch_seed_q_display(display_flag) &&
         gt_querymatch_seed_len_display(display_flag);
}

static int strcmp_ignore_ws(const char *s,const char *t)
{
  const char *sptr = s, *tptr = t;
  while (true)
  {
    if (isspace(*sptr))
    {
      sptr++;
    } else
    {
      if (isspace(*tptr))
      {
        tptr++;
      } else
      {
        if (*sptr < *tptr)
        {
          return -1;
        }
        if (*sptr > *tptr)
        {
          return 1;
        }
        if (*sptr == '\0')
        {
          return 0;
        }
        sptr++;
        tptr++;
      }
    }
  }
}

static const GtSEdisplayStruct *gt_display_arg_get(char *copyspace,
                                                   const char *str,
                                                   size_t cmplen)
{
  char *copy = NULL;
  const GtSEdisplayStruct *left = gt_display_arguments_table,
                          *right = gt_display_arguments_table +
                                   sizeof gt_display_arguments_table/
                                   sizeof gt_display_arguments_table[0] - 1;

  if (cmplen > 0)
  {
    copy = copyspace;
    memcpy(copy,str,sizeof *copy * cmplen);
    copy[cmplen] = '\0';
  }
  while (left <= right)
  {
    const GtSEdisplayStruct *mid = left + (right - left + 1)/2;
    const int cmp = strcmp_ignore_ws(copy == NULL ? str : copy,mid->name);

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

const unsigned int *gt_querymatch_display_order(GtUword *numcolumns,
                                                const GtSeedExtendDisplayFlag
                                                   *display_flag)
{
  gt_assert(display_flag != NULL);
  *numcolumns = display_flag->nextfree;
  return &display_flag->order[0];
}

static void gt_querymatch_display_flag_add(GtSeedExtendDisplayFlag
                                            *display_flag,
                                           GtSeedExtendDisplay_enum flag)
{
  const uint64_t mask = gt_display_mask(flag);

  gt_assert(display_flag != NULL);
  if ((display_flag->flags & mask) == 0)
  {
    display_flag->flags |= mask;
    gt_assert(flag <= GT_DISPLAY_LARGEST_FLAG);
    if (gt_display_arguments_table[gt_display_flag2index[flag]].incolumn)
    {
      gt_assert(display_flag->nextfree <= GT_DISPLAY_LARGEST_FLAG);
      display_flag->order[display_flag->nextfree++] = flag;
    }
  }
}

GtUword gt_querymatch_display_alignmentwidth(const GtSeedExtendDisplayFlag
                                                *display_flag)
{
  return (display_flag == NULL) ? 0 : display_flag->alignmentwidth;
}

GtUword gt_querymatch_trace_delta_display(const GtSeedExtendDisplayFlag
                                                *display_flag)
{
  return (display_flag == NULL) ? 0 : display_flag->trace_delta;
}

bool gt_querymatch_alignment_display(const GtSeedExtendDisplayFlag
                                     *display_flag)
{
  return (display_flag != NULL &&
          display_flag->alignmentwidth > 0) ? true : false;
}

#define GT_SE_ASSERT_DISPLAY_ID(PTR,ID)\
        if ((PTR) == NULL)\
        {\
          fprintf(stderr,"%s, %d: illegal identifier %s\n",\
                  __FILE__,__LINE__,ID);\
          exit(GT_EXIT_PROGRAMMING_ERROR);\
        }

static int gt_querymatch_display_flag_set(char *copyspace,
                                          GtWord *parameter,
                                          GtSeedExtendDisplayFlag *display_flag,
                                          const char *arg,
                                          GtError *err)
{
  const char *exclude_list[] = {"alignment","cigar",
                                "alignment","cigarX",
                                "cigar","cigarX",
                                "blast","custom",
                                "trace","alignment",
                                "trace","cigar",
                                "trace","cigarX",
                                "dtrace","alignment",
                                "dtrace","cigar",
                                "dtrace","cigarX",
                                "dtrace","trace",
                                "gfa2","blast",
                                "gfa2","alignment",
                                "gfa2","custom",
                                "gfa2","failed_seed",
                                "gfa2","seed_in_algn"};
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
  dstruct = gt_display_arg_get(copyspace,arg,cmplen);
  if (dstruct == NULL)
  {
    gt_error_set(err,"illegal identifier \"%s\" as argument of options "
                     "-outfmt, possible identifiers are: %s",arg,
                     GT_SE_POSSIBLE_DISPLAY_ARGS);
    return -1;
  }
  gt_querymatch_display_flag_add(display_flag,dstruct->flag);
  if (dstruct->flag == Gt_Seed_display)
  {
    gt_querymatch_display_flag_add(display_flag,Gt_Seed_len_display);
    gt_querymatch_display_flag_add(display_flag,Gt_Seed_s_display);
    gt_querymatch_display_flag_add(display_flag,Gt_Seed_q_display);
  }
  for (ex_idx = 0; ex_idx < numexcl; ex_idx+=2)
  {
    const GtSEdisplayStruct
      *dstruct0 = gt_display_arg_get(NULL,exclude_list[ex_idx],0),
      *dstruct1 = gt_display_arg_get(NULL,exclude_list[ex_idx+1],0);

    GT_SE_ASSERT_DISPLAY_ID(dstruct0,exclude_list[ex_idx]);
    GT_SE_ASSERT_DISPLAY_ID(dstruct1,exclude_list[ex_idx+1]);
    if ((display_flag->flags & gt_display_mask(dstruct0->flag)) &&
        (display_flag->flags & gt_display_mask(dstruct1->flag)))
    {
      gt_error_set(err,"argument \"%s\" and \"%s\" of option -outfmt exclude "
                       "each other",exclude_list[ex_idx],
                                    exclude_list[ex_idx+1]);
      return -1;
    }
  }
  return cmplen > 0 ? 1 : 0;
}

void gt_querymatch_display_flag_delete(GtSeedExtendDisplayFlag *display_flag)
{
  if (display_flag != NULL)
  {
    gt_free(display_flag);
  }
}

static void gt_querymatch_display_multi_flag_add(
                                  GtSeedExtendDisplayFlag *display_flag,
                                  const GtSeedExtendDisplay_enum *flag_enum,
                                  size_t numflags)
{
  size_t fidx;

  for (fidx = 0; fidx < numflags; fidx++)
  {
    gt_querymatch_display_flag_add(display_flag,flag_enum[fidx]);
  }
}

static bool gt_querymatch_display_args_contain(const GtStrArray *display_args,
                                               const char *keyword)
{
  GtUword idx;

  for (idx = 0; idx < gt_str_array_size(display_args); idx++)
  {
    if (strcmp_ignore_ws(gt_str_array_get(display_args,idx),keyword) == 0)
    {
      return true;
    }
  }
  return false;
}

static int gt_querymatch_options_order_check(
                 const GtSeedExtendDisplayFlag *display_flag,
                 size_t num_gfa2_default_flags,
                 GtError *err)
{
  GtUword idx, num_columns;
  bool trace_or_cigar_found = false;
  const unsigned int *order = gt_querymatch_display_order(&num_columns,
                                                          display_flag);

  if (!gt_querymatch_trace_display(display_flag) &&
      !gt_querymatch_cigar_display(display_flag))
  {
    gt_error_set(err,"for gfa2 output specify either trace or cigar as "
                     "argument of option -outfmt");
    return -1;
  }
  gt_assert((GtUword) num_gfa2_default_flags < num_columns);
  for (idx = num_gfa2_default_flags; idx < num_columns; idx++)
  {
    if (order[idx] == Gt_Trace_display || order[idx] == Gt_Cigar_display)
    {
      trace_or_cigar_found = true;
    } else
    {
      if (!trace_or_cigar_found)
      {
        const unsigned int display_arg_num = gt_display_flag2index[order[idx]];

        gt_error_set(err,"for gfa2 output in the list of argument to option "
                         "-outfmt the keyword \"%s\" must come after trace "
                         "or cigar",
                         gt_display_arguments_table[display_arg_num].name);
        return -1;
      }
    }
  }
  return 0;
}

GtSeedExtendDisplayFlag *gt_querymatch_display_flag_new(
                           const GtStrArray *display_args,
                           GtSeedExtendDisplaySetMode setmode,
                           GtError *err)
{
  GtSeedExtendDisplayFlag *display_flag = gt_malloc(sizeof *display_flag);
  char copyspace[GT_MAX_DISPLAY_FLAG_LENGTH+1];
  bool haserr = false;
  size_t num_gfa2_default_flags = 0;
  GtUword da_idx;
  const uint64_t trace_mask = gt_display_mask(Gt_Trace_display) |
                              gt_display_mask(Gt_Dtrace_display);
  /* required implications:
      mandatory: S_seqnum, Q_seqnum, S_start, Q_start
                 either editdist or score
      either S_len or S_end must be set
      either Q_len or Q_end must be set
      if Strand is not set, then order of Q_start and Q_end gives strand
  */
  display_flag->alignmentwidth = 0;
  display_flag->trace_delta = 0;
  display_flag->nextfree = 0;
  display_flag->flags = 0;
  if (setmode != GT_SEED_EXTEND_DISPLAY_SET_NO)
  {
    if (gt_querymatch_display_args_contain(display_args,"blast"))
    {
      GtSeedExtendDisplay_enum blast_flags[] =
      {
        Gt_Queryid_display,
        Gt_Subjectid_display,
        Gt_Identity_display,
        Gt_Alignmentlength_display,
        Gt_Mismatches_display,
        Gt_Gapopens_display,
        Gt_Q_start_display,
        Gt_Q_end_display,
        Gt_S_start_display,
        Gt_S_end_display,
        Gt_Evalue_display,
        Gt_Bitscore_display
      };
      gt_querymatch_display_multi_flag_add(display_flag,blast_flags,
                                           sizeof blast_flags/
                                           sizeof blast_flags[0]);
    } else
    {
      if (gt_querymatch_display_args_contain(display_args,"gfa2"))
      {
        GtSeedExtendDisplay_enum gfa2_flags[] =
        {
          Gt_S_seqnum_display,
          Gt_Q_seqnum_display,
          Gt_S_start_display,
          Gt_S_end_display,
          Gt_Q_start_display,
          Gt_Q_end_display
        };
        gt_querymatch_display_multi_flag_add(display_flag,gfa2_flags,
                                             sizeof gfa2_flags/
                                             sizeof gfa2_flags[0]);
        num_gfa2_default_flags = sizeof gfa2_flags/sizeof gfa2_flags[0];
      } else
      {
        if (!gt_querymatch_display_args_contain(display_args,"custom"))
        {
          if (setmode == GT_SEED_EXTEND_DISPLAY_SET_STANDARD)
          {
            GtSeedExtendDisplay_enum standard_flags[] =
            {
              Gt_S_len_display,
              Gt_S_seqnum_display,
              Gt_S_start_display,
              Gt_Strand_display,
              Gt_Q_len_display,
              Gt_Q_seqnum_display,
              Gt_Q_start_display,
              Gt_Score_display,
              Gt_Editdist_display,
              Gt_Identity_display
            };
            gt_querymatch_display_multi_flag_add(display_flag,standard_flags,
                                                 sizeof standard_flags/
                                                 sizeof standard_flags[0]);
          } else
          {

            gt_assert(setmode == GT_SEED_EXTEND_DISPLAY_SET_EXACT);
            GtSeedExtendDisplay_enum exact_flags[] =
            {
              Gt_S_len_display,
              Gt_S_seqnum_display,
              Gt_S_start_display,
              Gt_Strand_display,
              Gt_Q_len_display,
              Gt_Q_seqnum_display,
              Gt_Q_start_display
            };
            gt_querymatch_display_multi_flag_add(display_flag,exact_flags,
                                                 sizeof exact_flags/
                                                 sizeof exact_flags[0]);
          }
        }
      }
    }
  }
  for (da_idx = 0; da_idx < gt_str_array_size(display_args); da_idx++)
  {
    const char *da = gt_str_array_get(display_args,da_idx);
    GtWord parameter;
    int ret = gt_querymatch_display_flag_set(copyspace,&parameter,
                                             display_flag,da,err);

    if (ret == -1)
    {
      haserr = true;
      break;
    }
    if (ret == 1)
    {
      /* the only flag with a parameter is Gt_Alignment_display */
      if (parameter < 0)
      {
        const char *key;

        if (display_flag->flags & gt_display_mask(Gt_Alignment_display))
        {
          key = "alignment";
        } else
        {
          if (display_flag->flags & gt_display_mask(Gt_Trace_display))
          {
            key = "trace";
          } else
          {
            gt_assert(display_flag->flags & gt_display_mask(Gt_Dtrace_display));
            key = "dtrace";
          }
        }
        gt_error_set(err,"integer following \"%s=\" must be positive",key);
        haserr = true;
        break;
      }
      if (display_flag->flags & gt_display_mask(Gt_Alignment_display))
      {
        display_flag->alignmentwidth = (GtUword) parameter;
      } else
      {
        gt_assert (display_flag->flags & trace_mask);
        display_flag->trace_delta = (GtUword) parameter;
      }
    }
  }
  if (!haserr)
  {
    if ((display_flag->flags & gt_display_mask(Gt_Alignment_display)) &&
         display_flag->alignmentwidth == 0)
    {
      display_flag->alignmentwidth = GT_SEED_EXTEND_DEFAULT_ALIGNMENT_WIDTH;
    }
    if ((display_flag->flags & trace_mask) && display_flag->trace_delta == 0)
    {
      display_flag->trace_delta = GT_SEED_EXTEND_DEFAULT_TRACE_DELTA;
    }
  }
  if (!haserr && num_gfa2_default_flags > 0)
  {
    if (gt_querymatch_options_order_check(display_flag,
                                          num_gfa2_default_flags,err) != 0)
    {
      haserr = true;
    }
  }
  if (haserr)
  {
    gt_querymatch_display_flag_delete(display_flag);
    return NULL;
  }
  return display_flag;
}

void gt_querymatch_Options_output(FILE *stream,int argc,const char **argv,
                                  bool idhistout,GtUword minidentity,
                                  GtUword historysize)
{
  int idx;
  bool minid_out = false, history_out = false;

  fprintf(stream,"# Options:");
  for (idx = 1; idx < argc; idx++) {
    if (strcmp(argv[idx],"-minidentity") == 0) {
      minid_out = true;
    }
    if (strcmp(argv[idx],"-history") == 0) {
      history_out = true;
    }
    fprintf(stream," %s", argv[idx]);
  }
  if (idhistout)
  {
    if (!minid_out)
    {
      fprintf(stream," -minidentity " GT_WU,minidentity);
    }
    if (!history_out)
    {
      fprintf(stream," -history " GT_WU,historysize);
    }
  }
  fputc('\n',stream);
}

static void gt_querymatch_display_keyword_out(FILE *stream,const char *s)
{
  const char *sptr;

  for (sptr = s; *sptr != '\0'; sptr++)
  {
    if (*sptr == '.')
    {
      fprintf(stream,". ");
    } else
    {
      fputc(*sptr,stream);
    }
  }
}

void gt_querymatch_Fields_output(FILE *stream,
                                 const GtSeedExtendDisplayFlag *display_flag)
{
  const unsigned int *column_order;
  GtUword numcolumns, idx;

  gt_assert(display_flag != NULL);
  column_order = gt_querymatch_display_order(&numcolumns,display_flag);
  gt_assert(numcolumns > 0);
  fprintf(stream,"# Fields: ");
  gt_assert(numcolumns <= GT_DISPLAY_LARGEST_FLAG);
  for (idx = 0; idx < numcolumns; idx++)
  {
    unsigned int argnum, flag = column_order[idx];
    gt_assert(flag < sizeof gt_display_flag2index/
                     sizeof gt_display_flag2index[0]);
    argnum = gt_display_flag2index[flag];
    gt_assert(argnum < sizeof gt_display_arguments_table/
                       sizeof gt_display_arguments_table[0]);
    if (flag == Gt_Identity_display)
    {
      fprintf(stream,"%% %s",gt_display_arguments_table[argnum].name);
    } else
    {
      gt_querymatch_display_keyword_out(stream,
                                        gt_display_arguments_table[argnum].
                                          name);
    }
    fprintf(stream,"%s",idx < numcolumns - 1 ? ", " : "\n");
  }
}

const char *gt_querymatch_flag2name(GtSeedExtendDisplay_enum flag)
{
  gt_assert(flag <= GT_DISPLAY_LARGEST_FLAG);
  return gt_display_arguments_table[gt_display_flag2index[flag]].name;
}

GtStrArray *gt_querymatch_read_Fields_line(const char *line_ptr)
{
  const char *header = "# Fields:", *ptr, *last_start;
  const size_t header_len = strlen(header);
  GtStrArray *fields;

  if (strncmp(header,line_ptr,header_len) != 0)
  {
    return NULL;
  }
  fields = gt_str_array_new();
  for (last_start = ptr = line_ptr + header_len + 1; /* Nothing */; ptr++)
  {
    if (*ptr == ',' || *ptr == '\0')
    {
      GtUword len_arg;

      if (*last_start == '%')
      {
        last_start += 2;
      }
      len_arg = (size_t) (ptr - last_start);
      gt_str_array_add_cstr_nt(fields, last_start, len_arg);
      last_start = ptr + 2;
    }
    if (*ptr == '\0')
    {
      break;
    }
  }
  return fields;
}
