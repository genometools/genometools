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

#include "core/ma_api.h"
#include "core/str_api.h"
#include "core/encseq.h"
#include "match/querymatch.h"
#include "match/seed-extend.h"
#include "match/seed-extend-iter.h"

struct GtSeedextendMatchIterator
{
  GtStr *ii, *qii;
  bool mirror, bias_parameters, seqlength_display;
  GtEncseq *aencseq, *bencseq;
  GtUword errorpercentage, history_size;
  const char *matchfilename;
  GtStr *line_buffer;
  uint64_t linenum;
  bool has_seedline;
  GtUword seedpos1,
          seedpos2,
          seedlen;
  FILE *inputfileptr;
  GtUword currentmatchindex;
  GtQuerymatch *currentmatch, *querymatchptr;
  GtQuerymatchoutoptions *querymatchoutoptions;
  GtArrayGtQuerymatch querymatch_table;
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
  FILE *options_line_inputfileptr;
  int had_err = 0;
  GtSeedextendMatchIterator *semi = gt_malloc(sizeof *semi);
  GtEncseqLoader *encseq_loader = NULL;

  semi->ii = gt_str_new(),
  semi->qii = gt_str_new();
  semi->mirror = false;
  semi->bias_parameters = false;
  semi->seqlength_display = false;
  semi->aencseq = semi->bencseq = NULL;
  semi->errorpercentage = 0;
  semi->history_size = 0;
  semi->matchfilename = gt_str_get(matchfilename);
  semi->has_seedline = false;
  semi->line_buffer = NULL;
  semi->inputfileptr = NULL;
  semi->querymatchptr = gt_querymatch_new();
  semi->currentmatchindex = GT_UWORD_MAX;
  semi->currentmatch = NULL;
  semi->querymatchoutoptions = NULL;
  semi->seedpos1 = semi->seedpos2 = semi->seedlen = GT_UWORD_MAX;
  GT_INITARRAY(&semi->querymatch_table,GtQuerymatch);
  options_line_inputfileptr = fopen(semi->matchfilename, "r");
  if (options_line_inputfileptr == NULL)
  {
    gt_error_set(err, "file %s does not exist", semi->matchfilename);
    had_err = -1;
  }
  if (!had_err)
  {
    GtStr *options_line_buffer = gt_str_new();
    /* read first line and evaluate tokens */
    if (gt_str_read_next_line(options_line_buffer,
                              options_line_inputfileptr) != EOF)
    {
      char *tok, *lineptr = gt_str_get(options_line_buffer);
      bool parse_ii = false, parse_qii = false, parse_minid = false,
           parse_history = false;;

      while (!had_err && (tok = strsep(&lineptr," ")) != NULL)
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
          continue;
        }
        if (strcmp(tok, "-qii") == 0)
        {
          parse_qii = true;
          continue;
        }
        if (strcmp(tok, "-minidentity") == 0)
        {
          parse_minid = true;
          continue;
        }
        if (strcmp(tok, "-history") == 0)
        {
          parse_history = true;
          continue;
        }
        if (strcmp(tok, "-mirror") == 0)
        {
          semi->mirror = true; /* found -mirror option */
        }
        if (strcmp(tok, "-bias-parameters") == 0)
        {
          semi->bias_parameters = true;
        }
        if (strcmp(tok, "-seqlength-display") == 0)
        {
          semi->seqlength_display = true;
        }
      }
      if (!had_err)
      {
        if (gt_str_length(semi->ii) == 0UL)
        {
          gt_error_set(err, "need output of option string "
                            "(run gt seed_extend with -v or -verify)");
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
    }
    gt_str_delete(options_line_buffer);
  }
  if (options_line_inputfileptr != NULL)
  {
    fclose(options_line_inputfileptr);
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
    /* ignore comment lines; but print seeds if -seed-display is set */
    gt_assert(line_ptr != NULL);
    if (line_ptr[0] != '\n')
    {
      if (line_ptr[0] == '#')
      {
        const char *seedwordptr = strstr(line_ptr, "seed:");

        if (seedwordptr != NULL)
        {
          int num = sscanf(seedwordptr + sizeof ("seed:"),
                           GT_WU " " GT_WU " " GT_WU,
                           &semi->seedpos1,&semi->seedpos2,
                           &semi->seedlen);
          semi->has_seedline = num == 3 ? true : false;
        } else
        {
          semi->has_seedline = false;
          semi->seedpos1 = semi->seedpos2 = semi->seedlen = GT_UWORD_MAX;
        }
      } else
      {
        if (gt_querymatch_read_line(semi->querymatchptr,
                                    semi->seqlength_display,
                                    line_ptr,
                                    selfmatch,
                                    semi->seedpos1,
                                    semi->seedpos2,
                                    semi->seedlen,
                                    semi->aencseq,
                                    semi->bencseq))
        {
          gt_str_reset(semi->line_buffer);
          semi->seedpos1 = semi->seedpos2 = semi->seedlen = GT_UWORD_MAX;
          return semi->querymatchptr;
        }
      }
      gt_str_reset(semi->line_buffer);
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

bool gt_seedextend_match_iterator_seqlength_display(
                        const GtSeedextendMatchIterator *semi)
{
  gt_assert(semi != NULL);
  return semi->seqlength_display;
}

bool gt_seedextend_match_iterator_has_seedline(
                        const GtSeedextendMatchIterator *semi)
{
  gt_assert(semi != NULL);
  if (semi->currentmatch != NULL)
  {
    return gt_querymatch_has_seed(semi->currentmatch);
  }
  return semi->has_seedline;
}

GtQuerymatch *gt_seedextend_match_iterator_querymatch_ptr(
                        GtSeedextendMatchIterator *semi)
{
  gt_assert(semi != NULL);
  return semi->querymatchptr;
}

GtUword gt_seedextend_match_iterator_seedlen(
                          const GtSeedextendMatchIterator *semi)
{
  gt_assert(semi != NULL && semi->has_seedline);
  return semi->seedlen;
}

GtUword gt_seedextend_match_iterator_seedpos1(
                         const GtSeedextendMatchIterator *semi)
{
  gt_assert(semi != NULL && semi->has_seedline);
  return semi->seedpos1;
}

GtUword gt_seedextend_match_iterator_seedpos2(
                        const GtSeedextendMatchIterator *semi)
{
  gt_assert(semi != NULL && semi->has_seedline);
  return semi->seedpos2;
}

void gt_seedextend_match_iterator_display_set(GtSeedextendMatchIterator *semi,
                                              unsigned int display_flag)
{
  gt_assert(semi != NULL);
  gt_querymatch_display_set(semi->querymatchptr,display_flag);
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

void gt_seedextend_match_iterator_querymatchoutoptions_set(
                    GtSeedextendMatchIterator *semi,
                    bool generatealignment,
                    bool showeoplist,
                    GtUword alignmentwidth,
                    bool always_polished_ends,
                    unsigned int display_flag)
{
  double matchscore_bias = GT_DEFAULT_MATCHSCORE_BIAS;

  semi->querymatchoutoptions
    = gt_querymatchoutoptions_new(generatealignment,showeoplist,alignmentwidth);
  if (gt_seedextend_match_iterator_bias_parameters(semi))
  {
    matchscore_bias = gt_greedy_dna_sequence_bias_get(semi->aencseq);
  }
  gt_querymatchoutoptions_for_align_only(semi->querymatchoutoptions,
                            semi->errorpercentage,
                            matchscore_bias,
                            gt_seedextend_match_iterator_history_size(semi),
                            always_polished_ends,
                            display_flag);
  gt_querymatch_outoptions_set(semi->querymatchptr,semi->querymatchoutoptions);
}
