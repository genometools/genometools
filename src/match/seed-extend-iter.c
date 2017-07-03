/*
  Copyright (c) 2015 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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

#include <float.h>
#include "core/ma_api.h"
#include "core/str_api.h"
#include "core/encseq.h"
#include "match/querymatch.h"
#include "match/seed-extend.h"
#include "match/seed-extend-iter.h"

struct GtSeedextendMatchIterator
{
  GtStr *ii, *qii;
  bool mirror, bias_parameters;
  GtEncseq *aencseq, *bencseq;
  GtUword errorpercentage, history_size;
  const char *matchfilename;
  GtStr *line_buffer;
  uint64_t linenum;
  FILE *inputfileptr;
  double  evalue, bitscore;
  GtUword currentmatchindex;
  GtQuerymatch *currentmatch, *querymatchptr;
  GtQuerymatchoutoptions *querymatchoutoptions;
  GtArrayGtQuerymatch querymatch_table;
  GtSeedExtendDisplayFlag *in_display_flag;
  GtStr *saved_options_line;
  GtUword trace_delta;
  bool missing_fields_line;
};

void gt_seedextend_match_iterator_delete(GtSeedextendMatchIterator *semi)
{
  if (semi == NULL)
  {
    return;
  }
  gt_str_delete(semi->ii);
  gt_str_delete(semi->qii);
  gt_encseq_delete(semi->aencseq);
  gt_encseq_delete(semi->bencseq);
  gt_str_delete(semi->line_buffer);
  gt_querymatch_delete(semi->querymatchptr);
  GT_FREEARRAY(&semi->querymatch_table,GtQuerymatch);
  gt_querymatchoutoptions_delete(semi->querymatchoutoptions);
  if (semi->inputfileptr != NULL)
  {
    fclose(semi->inputfileptr);
  }
  gt_querymatch_display_flag_delete(semi->in_display_flag);
  gt_str_delete(semi->saved_options_line);
  gt_free(semi);
}

#if defined(_WIN32) || defined(_WIN64)
char *strsep(char** stringp, const char *delim)
{
  char *p, *start = *stringp;

  p = (start != NULL) ? strpbrk(start, delim) : NULL;
  if (p == NULL)
  {
    *stringp = NULL;
  } else
  {
    *p = '\0';
    *stringp = p + 1;
  }
  return start;
}
#endif

/* Parse encseq input indices from -ii and -qii options in 1st line of file. */
GtSeedextendMatchIterator *gt_seedextend_match_iterator_new(
                                            const GtStr *matchfilename,
                                            GtError *err)
{
  FILE *defline_infp;
  int had_err = 0;
  GtSeedextendMatchIterator *semi = gt_malloc(sizeof *semi);
  GtEncseqLoader *encseq_loader = NULL;

  semi->ii = gt_str_new(),
  semi->qii = gt_str_new();
  semi->mirror = false;
  semi->bias_parameters = false;
  semi->aencseq = semi->bencseq = NULL;
  semi->errorpercentage = 0;
  semi->history_size = 0;
  semi->matchfilename = gt_str_get(matchfilename);
  semi->line_buffer = NULL;
  semi->inputfileptr = NULL;
  semi->querymatchptr = gt_querymatch_new();
  semi->currentmatchindex = GT_UWORD_MAX;
  semi->evalue = DBL_MAX;
  semi->bitscore = DBL_MAX;
  semi->currentmatch = NULL;
  semi->querymatchoutoptions = NULL;
  semi->in_display_flag = NULL;
  semi->trace_delta = GT_SEED_EXTEND_DEFAULT_TRACE_DELTA;
  semi->saved_options_line = NULL;
  GT_INITARRAY(&semi->querymatch_table,GtQuerymatch);
  defline_infp = fopen(semi->matchfilename, "r");
  if (defline_infp == NULL)
  {
    gt_error_set(err, "file %s does not exist", semi->matchfilename);
    had_err = -1;
  }
  if (!had_err)
  {
    GtStr *options_line_buffer = gt_str_new();
    /* read first line and evaluate tokens */
    if (gt_str_read_next_line(options_line_buffer,defline_infp) != EOF)
    {
      char *tok, *line_ptr = gt_str_get(options_line_buffer);
      bool parse_ii = false, parse_qii = false, parse_minid = false,
           parse_history = false, parse_outfmt = false;

      semi->saved_options_line = gt_str_clone(options_line_buffer);
      while (!had_err && (tok = strsep(&line_ptr," ")) != NULL)
      {
        if (parse_minid || parse_history)
        {
          GtWord readgtword;

          if (sscanf(tok,GT_WD,&readgtword) != 1 ||
              readgtword < 0  || (parse_minid && readgtword > 99) ||
              (parse_history && readgtword > 64))
          {
            gt_error_set(err,"cannot parse argument for option -%s from "
                             "first line of file %s",
                             parse_minid ? "minidentity"
                                         : "history",
                             semi->matchfilename);
            had_err = -1;
          }
          if (parse_minid)
          {
            semi->errorpercentage = 100 - readgtword;
            parse_minid = false;
          } else
          {
            semi->history_size = (GtUword) readgtword;
            parse_history = false;
          }
          continue;
        }
        if (parse_outfmt)
        {
          GtWord this_trace_delta;

          if (sscanf(tok,"trace=" GT_WD,&this_trace_delta) == 1 ||
              sscanf(tok,"dtrace=" GT_WD,&this_trace_delta) == 1)
          {
            if (this_trace_delta < 0)
            {
              gt_error_set(err,"value of trace= cannot be negative");
              had_err = -1;
            }
            semi->trace_delta = (GtUword) this_trace_delta;
          }
        }
        if (parse_ii)
        {
          gt_str_set(semi->ii, tok);
          parse_ii = false;
          continue;
        }
        if (parse_qii)
        {
          gt_str_set(semi->qii, tok);
          parse_qii = false;
          continue;
        }
        if (strcmp(tok, "-ii") == 0)
        {
          parse_ii = true;
          parse_outfmt = false;
          continue;
        }
        if (strcmp(tok, "-qii") == 0)
        {
          parse_qii = true;
          parse_outfmt = false;
          continue;
        }
        if (strcmp(tok, "-minidentity") == 0)
        {
          parse_minid = true;
          parse_outfmt = false;
          continue;
        }
        if (strcmp(tok, "-history") == 0)
        {
          parse_history = true;
          parse_outfmt = false;
          continue;
        }
        if (strcmp(tok, "-mirror") == 0)
        {
          parse_outfmt = false;
          semi->mirror = true; /* found -mirror option */
        }
        if (strcmp(tok, "-bias-parameters") == 0)
        {
          parse_outfmt = false;
          semi->bias_parameters = true;
        }
        if (strcmp(tok, "-outfmt") == 0)
        {
          parse_outfmt = true;
        }
      }
      if (!had_err)
      {
        if (gt_str_length(semi->ii) == 0UL)
        {
          gt_error_set(err, "missing option string");
          had_err = -1;
        }
        if (semi->errorpercentage == 0)
        {
          gt_error_set(err,"missing option -minidentity in first line of file "
                           "%s",semi->matchfilename);
          had_err = -1;
        }
      }
    } else
    {
      gt_error_set(err, "file %s is empty", semi->matchfilename);
      had_err = -1;
    };
    gt_str_delete(options_line_buffer);
  }
  if (!had_err)
  {
    GtStr *fieldsline_buffer = gt_str_new();
    GtStrArray *fields = NULL;

    while (gt_str_read_next_line(fieldsline_buffer,defline_infp) != EOF &&
           fields == NULL)
    {
      char *line_ptr = gt_str_get(fieldsline_buffer);
      fields = gt_querymatch_read_Fields_line(line_ptr);

      gt_str_reset(fieldsline_buffer);
    }
    gt_str_delete(fieldsline_buffer);
    if (fields != NULL)
    {
      const GtSeedExtendDisplaySetMode setmode = GT_SEED_EXTEND_DISPLAY_SET_NO;

      semi->in_display_flag = gt_querymatch_display_flag_new(fields,setmode,
                                                             err);
      if (semi->in_display_flag == NULL)
      {
        had_err = -1;
      }
      gt_str_array_delete(fields);
    }
  }
  if (defline_infp != NULL)
  {
    fclose(defline_infp);
  }
  /* Load encseqs */
  if (!had_err)
  {
    encseq_loader = gt_encseq_loader_new();
    gt_encseq_loader_enable_autosupport(encseq_loader);
    semi->aencseq
      = gt_encseq_loader_load(encseq_loader, gt_str_get(semi->ii), err);
    if (semi->aencseq == NULL)
    {
      had_err = -1;
    }
  }
  if (!had_err)
  {
    if (gt_str_length(semi->qii) != 0)
    {
      semi->bencseq
        = gt_encseq_loader_load(encseq_loader, gt_str_get(semi->qii), err);
    } else {
      semi->bencseq = gt_encseq_ref(semi->aencseq);
    }
    if (semi->bencseq == NULL)
    {
      had_err = -1;
      gt_encseq_delete(semi->aencseq);
    }
  }
  gt_encseq_loader_delete(encseq_loader);
  if (!had_err)
  {
    gt_assert(semi->bencseq != NULL);
    semi->line_buffer = gt_str_new();
    semi->linenum = 0;
    semi->inputfileptr = fopen(semi->matchfilename, "rb");
    if (semi->inputfileptr == NULL)
    {
      gt_error_set(err, "file %s does not exist", semi->matchfilename);
      had_err = true;
    }
  }
  if (had_err)
  {
    gt_seedextend_match_iterator_delete(semi);
    return NULL;
  }
  return semi;
}

GtQuerymatch *gt_seedextend_match_iterator_next(GtSeedextendMatchIterator *semi)
{
  bool selfmatch;

  gt_assert(semi != NULL);
  if (semi->currentmatchindex < GT_UWORD_MAX)
  {
    if (semi->currentmatchindex == semi->querymatch_table.nextfreeGtQuerymatch)
    {
      semi->currentmatch = NULL;
    } else
    {
      semi->currentmatch = gt_querymatch_table_get(&semi->querymatch_table,
                                                   semi->currentmatchindex++);
    }
    return semi->currentmatch;
  }
  selfmatch = semi->aencseq == semi->bencseq ? true : false;
  while (true)
  {
    const char *line_ptr;

    semi->linenum++;
    if (gt_str_read_next_line(semi->line_buffer, semi->inputfileptr) == EOF)
    {
      break;
    }
    line_ptr = gt_str_get(semi->line_buffer);
    gt_assert(line_ptr != NULL);
    if (line_ptr[0] != '\n' && line_ptr[0] != '#')
    {
      if (semi->in_display_flag == NULL)
      {
        return NULL;
      }
      gt_querymatch_read_line(semi->querymatchptr,
                              &semi->evalue,
                              &semi->bitscore,
                              line_ptr,
                              semi->in_display_flag,
                              selfmatch,
                              semi->aencseq,
                              semi->bencseq);
      gt_str_reset(semi->line_buffer);
      return semi->querymatchptr;
    }
    gt_str_reset(semi->line_buffer);
  }
  return NULL;
}

const GtEncseq *gt_seedextend_match_iterator_aencseq(
                        const GtSeedextendMatchIterator *semi)
{
  gt_assert(semi != NULL);
  return semi->aencseq;
}

const GtEncseq *gt_seedextend_match_iterator_bencseq(
                        const GtSeedextendMatchIterator *semi)
{
  gt_assert(semi != NULL);
  return semi->bencseq;
}

GtUword gt_seedextend_match_iterator_history_size(
                        const GtSeedextendMatchIterator *semi)
{
  gt_assert(semi != NULL);
  return semi->history_size;
}

GtUword gt_seedextend_match_iterator_errorpercentage(
                        const GtSeedextendMatchIterator *semi)
{
  gt_assert(semi != NULL);
  return semi->errorpercentage;
}

bool gt_seedextend_match_iterator_bias_parameters(
                        const GtSeedextendMatchIterator *semi)
{
  gt_assert(semi != NULL);
  return semi->bias_parameters;
}

bool gt_seedextend_match_iterator_has_seed(
                        const GtSeedextendMatchIterator *semi)
{
  gt_assert(semi != NULL);
  return gt_querymatch_has_seed(semi->in_display_flag);
}

bool gt_seedextend_match_iterator_has_cigar(
                        const GtSeedextendMatchIterator *semi)
{
  gt_assert(semi != NULL);
  return gt_querymatch_cigar_display(semi->in_display_flag) ||
         gt_querymatch_cigarX_display(semi->in_display_flag);
}

GtUword gt_seedextend_match_iterator_trace_delta(
                        const GtSeedextendMatchIterator *semi)
{
  gt_assert(semi != NULL);
  return (gt_querymatch_trace_display(semi->in_display_flag) ||
          gt_querymatch_dtrace_display(semi->in_display_flag))
           ? semi->trace_delta : 0;
}

bool gt_seedextend_match_iterator_dtrace(const GtSeedextendMatchIterator *semi)
{
  return gt_querymatch_dtrace_display(semi->in_display_flag) ? true : false;
}

double gt_seedextend_match_iterator_evalue(const GtSeedextendMatchIterator
                                              *semi)
{
  gt_assert(semi != NULL);
  return semi->evalue;
}

double gt_seedextend_match_iterator_bitscore(
                  const GtSeedextendMatchIterator *semi)
{
  gt_assert(semi != NULL);
  return semi->bitscore;
}

void gt_seedextend_match_iterator_verify_alignment_set(
                                  GtSeedextendMatchIterator *semi)
{
  gt_querymatch_verify_alignment_set(semi->querymatchptr);
}

GtUword gt_seedextend_match_iterator_all_sorted(GtSeedextendMatchIterator *semi,
                                                bool ascending)

{
  GtQuerymatch *querymatchptr;
  gt_assert(semi != NULL);

  while ((querymatchptr = gt_seedextend_match_iterator_next(semi)) != NULL)
  {
    gt_querymatch_table_add(&semi->querymatch_table,querymatchptr);
  }
  gt_querymatch_table_sort(&semi->querymatch_table,ascending);
  semi->currentmatchindex = 0;
  semi->currentmatch = NULL;
  return semi->querymatch_table.nextfreeGtQuerymatch;
}

GtQuerymatch *gt_seedextend_match_iterator_get(
                            const GtSeedextendMatchIterator *semi,
                            GtUword idx)
{
  gt_assert(semi != NULL && idx < semi->querymatch_table.nextfreeGtQuerymatch);

  return gt_querymatch_table_get(&semi->querymatch_table,idx);
}

int gt_seedextend_match_iterator_querymatchoutoptions_set(
                    GtSeedextendMatchIterator *semi,
                    bool always_polished_ends,
                    GtExtendCharAccess a_extend_char_access,
                    GtExtendCharAccess b_extend_char_access,
                    const GtSeedExtendDisplayFlag *out_display_flag,
                    GtError *err)
{
  double matchscore_bias = GT_DEFAULT_MATCHSCORE_BIAS;

  semi->querymatchoutoptions
    = gt_querymatchoutoptions_new(out_display_flag,gt_str_get(semi->ii),err);
  if (semi->querymatchoutoptions == NULL)
  {
    return -1;
  }
  gt_assert(semi->in_display_flag != NULL && out_display_flag != NULL);
  if (gt_querymatch_cigar_display(semi->in_display_flag) &&
      gt_querymatch_cigarX_display(out_display_flag))
  {
    gt_error_set(err,"match file with alignments in cigar format cannot be "
                     "converted to cigarX format");
    return -1;
  }
  if (gt_seedextend_match_iterator_bias_parameters(semi))
  {
    matchscore_bias = gt_greedy_dna_sequence_bias_get(semi->aencseq);
  }
  gt_querymatchoutoptions_for_align_only(semi->querymatchoutoptions,
                            semi->errorpercentage,
                            matchscore_bias,
                            gt_seedextend_match_iterator_history_size(semi),
                            always_polished_ends,
                            a_extend_char_access,
                            b_extend_char_access,
                            out_display_flag);
  gt_querymatch_outoptions_set(semi->querymatchptr,semi->querymatchoutoptions);
  return 0;
}

const char *gt_seedextend_match_iterator_Options_line(
                    const GtSeedextendMatchIterator *semi)
{
  gt_assert(semi != NULL && semi->saved_options_line != NULL);
  return gt_str_get(semi->saved_options_line);
}
