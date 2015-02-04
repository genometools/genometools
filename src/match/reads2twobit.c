/*
  Copyright (c) 2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#include <ctype.h>
#include <errno.h>
#include <limits.h>

#include "core/array.h"
#include "core/basename_api.h"
#include "core/encseq.h"
#include "core/fa.h"
#include "core/fileutils.h"
#include "core/filelengthvalues.h"
#include "core/intbits.h"
#include "core/log_api.h"
#include "core/minmax.h"
#include "core/parseutils_api.h"
#include "core/qsort_r_api.h"
#include "core/splitter_api.h"
#include "core/desc_buffer.h"
#include "core/xansi_api.h"
#include "core/undef_api.h"
#include "match/hplstore.h"
#include "match/reads_libraries_table.h"
#include "match/reads2twobit.h"

#define GT_READS2TWOBIT_ALPHASIZE 4U

typedef struct {
  bool paired;
  GtStr *filename1;
  GtStr *filename2;
  GtUword insertlength;
  GtUword stdev; /* 0 if unknown */
  GtUword total_filelength;
  GtUword total_seqlength;
  GtUword first_seqnum;
  GtUword nofseqs;
} GtReadsLibraryInfo;

struct GtReads2Twobit
{
  GtStr *indexname;
  GtArray *collection;
  GtTwobitencoding *twobitencoding;
  GtUword chardistri[GT_READS2TWOBIT_ALPHASIZE];
  GtUword *seppos;
  GtUword nofseqs;
  GtUword seqlen_eqlen, seqlen_max, seqlen_min;
  GtUword total_seqlength;
  GtTwobitencoding current_sepcode;
  GtUword invalid_sequences, invalid_total_length;
  char phredbase, lowqual;
  GtUword maxlow;
  bool has_paired;
  bool use_rle;
  GtDescBuffer *descs;
  FILE *descsfp;
  bool clipdes;
  GtUword longestdesc;
  GtUword n_descs;
  GtHplstore *hplengths;
  double approx_avhlen;
};

GtReads2Twobit* gt_reads2twobit_new(GtStr *indexname)
{
  GtReads2Twobit *r2t;
  gt_assert(indexname != NULL);
  r2t = gt_malloc(sizeof (*r2t));
  r2t->indexname = indexname;
  r2t->collection = gt_array_new(sizeof (GtReadsLibraryInfo));
  r2t->twobitencoding = NULL;
  r2t->seppos = NULL;
  r2t->invalid_sequences = 0;
  r2t->invalid_total_length = 0;
  r2t->seqlen_eqlen = 0;
  r2t->seqlen_max = 0;
  r2t->seqlen_min = 0;
  r2t->total_seqlength = 0;
  r2t->nofseqs = 0;
  r2t->phredbase = (char)33;
  r2t->maxlow = GT_UNDEF_UWORD;
  r2t->lowqual = 0;
  r2t->has_paired = false;
  r2t->use_rle = false;
  r2t->hplengths = NULL;
  r2t->approx_avhlen = 0.0;
  r2t->descs = NULL;
  r2t->descsfp = NULL;
  r2t->clipdes = true;
  r2t->longestdesc = 0;
  r2t->n_descs = 0;
  return r2t;
}

void gt_reads2twobit_use_rle(GtReads2Twobit *r2t)
{
  gt_assert(r2t != NULL);
  r2t->use_rle = true;
}

void gt_reads2twobit_delete(GtReads2Twobit *r2t)
{
  if (r2t != NULL)
  {
    GtUword lnum, noflibs;
    noflibs = gt_array_size(r2t->collection);
    for (lnum = 0; lnum < noflibs; lnum++)
    {
      GtReadsLibraryInfo *rli;
      rli = gt_array_get(r2t->collection, lnum);
      gt_str_delete(rli->filename1);
      if (rli->filename2 != NULL)
        gt_str_delete(rli->filename2);
    }
    gt_array_delete(r2t->collection);
    gt_desc_buffer_delete(r2t->descs);
    gt_fa_xfclose(r2t->descsfp);
    gt_hplstore_delete(r2t->hplengths);
    gt_free(r2t->twobitencoding);
    gt_free(r2t->seppos);
    gt_free(r2t);
  }
}

bool gt_reads2twobit_has_paired(GtReads2Twobit *r2t)
{
  gt_assert(r2t != NULL);
  return r2t->has_paired;
}

void gt_reads2twobit_add(GtReads2Twobit *r2t, bool paired,
    const GtStr *filename1, const GtStr *filename2, GtUword insertlength,
    GtUword stdev)
{
  GtReadsLibraryInfo rli;
  gt_assert(r2t != NULL);
  gt_assert(filename1 != NULL);
  rli.paired = paired;
  rli.filename1 = gt_str_clone(filename1);
  if (filename2 != NULL)
    rli.filename2 = gt_str_clone(filename2);
  else
    rli.filename2 = NULL;
  rli.insertlength = insertlength;
  rli.stdev = stdev;
  rli.total_seqlength = 0;
  rli.first_seqnum = 0;
  rli.nofseqs = 0;
  rli.total_filelength = (GtUword)gt_file_size(gt_str_get(filename1));
  if (filename2 != NULL)
    rli.total_filelength += (GtUword)gt_file_size(gt_str_get(filename2));
  gt_array_add(r2t->collection, rli);
}

int gt_reads2twobit_add_library(GtReads2Twobit *r2t, const GtStr *libspec,
    GtError *err)
{
  GtSplitter *s1, *s2;
  char *libspec_copy;
  int had_err = 0;

  gt_assert(r2t != NULL);
  gt_assert(libspec != NULL);
  gt_assert(gt_str_length(libspec) > 0);
  libspec_copy = gt_malloc(sizeof (*libspec_copy) *
      (gt_str_length(libspec) + 1UL));
  (void)strcpy(libspec_copy, gt_str_get(libspec));
  s1 = gt_splitter_new();
  gt_splitter_split(s1, libspec_copy, gt_str_length(libspec),
      GT_READS2TWOBIT_LIBSPECSEP);
  gt_log_log("reads2twobit add library: %s", gt_str_get(libspec));
  if (gt_splitter_size(s1) == 1UL)
  {
    gt_reads2twobit_add(r2t, false, libspec, NULL, 0, 0);
  }
  else if (gt_splitter_size(s1) < 4UL)
  {
    GtStr *filename1, *filename2;
    GtUword insertlength = 0, stdev = 0;
    char *insertspec, *insertspec_copy;
    r2t->has_paired = true;
    filename1 = gt_str_new_cstr(gt_splitter_get_token(s1, 0));
    if (gt_splitter_size(s1) == 3UL)
    {
      filename2 = gt_str_new_cstr(gt_splitter_get_token(s1, 1UL));
      insertspec = gt_splitter_get_token(s1, 2UL);
    }
    else
    {
      filename2 = NULL;
      insertspec = gt_splitter_get_token(s1, 1UL);
    }
    s2 = gt_splitter_new();
    insertspec_copy = gt_malloc(sizeof (*insertspec_copy) *
        (strlen(insertspec) + 1UL));
    (void)strcpy(insertspec_copy, insertspec);
    gt_splitter_split(s2, insertspec_copy, (GtUword)strlen(insertspec),
      GT_READS2TWOBIT_INSERTSEP);
    if (gt_splitter_size(s2) <= 2UL)
    {
      had_err = gt_parse_uword(&insertlength, gt_splitter_get_token(s2, 0));
      if (gt_splitter_size(s2) == 2UL && !had_err)
      {
        had_err = gt_parse_uword(&stdev, gt_splitter_get_token(s2, 1UL));
      }
    }
    if (gt_splitter_size(s2) > 2UL || had_err)
    {
      gt_error_set(err, "insert specification not valid: %s\nthe correct "
          "syntax is \"insertlength[%cstdev]\"", insertspec,
          GT_READS2TWOBIT_INSERTSEP);
    }
    if (!had_err)
    {
      gt_reads2twobit_add(r2t, true, filename1, filename2, insertlength, stdev);
    }
    gt_str_delete(filename1);
    gt_str_delete(filename2);
    gt_free(insertspec_copy);
    gt_splitter_delete(s2);
  }
  else
  {
    gt_error_set(err, "library specification not valid: %s\nthe correct "
        "syntax is \"filename[[%cfilename2]%cinsertlength[%cstdev]]\"\n"
        "(filenames are not allowed to contain \"%c\"",
        gt_str_get(libspec),
      GT_READS2TWOBIT_LIBSPECSEP, GT_READS2TWOBIT_LIBSPECSEP,
      GT_READS2TWOBIT_INSERTSEP, GT_READS2TWOBIT_LIBSPECSEP);
    had_err = -1;
  }
  gt_free(libspec_copy);
  gt_splitter_delete(s1);
  return had_err;
}

void gt_reads2twobit_use_phred64(GtReads2Twobit *r2t)
{
  gt_assert(r2t != NULL);
  r2t->phredbase = (char)64;
}

void gt_reads2twobit_set_quality_filter(GtReads2Twobit *r2t,
    GtUword maxlow, char lowqual)
{
  r2t->maxlow = maxlow;
  r2t->lowqual = lowqual;
}

#define GT_READS2TWOBIT_CODE_UNDEF ((GtTwobitencoding)ULONG_MAX)

typedef struct {
  GtTwobitencoding *tbe_next;
  GtTwobitencoding kmercode;
  unsigned short codepos;
  GtUword chardistri[GT_READS2TWOBIT_ALPHASIZE];
  GtUword globalpos;
  GtUword nofseqs;
  GtUword seppos_nextfree;
  GtUword seqlen_max;
  GtUword seqlen_min;
  GtUword seqlen_first;
} GtReads2TwobitEncodeInfo;

#define GT_READS2TWOBIT_INIT_ENCODE_INFO(EI, TBE)\
  memset(&(EI), 0, sizeof (EI));\
  (EI).tbe_next = (TBE)

#define GT_READS2TWOBIT_COPY_ENCODE_INFO(EISRC, EIDEST)\
  memcpy(&(EIDEST), &(EISRC), sizeof (EISRC))

typedef struct {
  GtTwobitencoding char2code[UCHAR_MAX + 1];
  GtUword inputfiles_totallength;
  GtUword seqlen, seqlen_mate, exp_qlen;
  GtReads2TwobitEncodeInfo current, backup;
  bool varlen_mode, invalid_mode;
  GtUword invalid_sequences, invalid_total_length;
  GtUword *seppos;
  GtUword seppos_alloc;
  char phredbase, *qbuf, *qbuf2, lowqual;
  size_t qbuf_size, qbuf2_size;
  GtStr *dbuf, *dbuf2;
  GtUword maxlow;
  bool use_rle;
  GtTwobitencoding prevcode;
  GtHplstore *hplengths;
  uint8_t hplength;
  GtUword hsum, nofh;
  GtDescBuffer *descs;
  FILE *descsfp;
  bool clipdes;
  GtUword *longestdesc;
  GtUword *n_descs;
} GtReads2TwobitEncodeState;

static int gt_reads2twobit_rli_cmp(const void *a, const void *b)
{
  const GtReadsLibraryInfo *rli_a = a, *rli_b = b;
  return (int)(rli_b->paired - rli_a->paired);
}

static void gt_reads2twobit_init_encode(GtReads2Twobit *r2t,
    GtReads2TwobitEncodeState *state)
{
  const GtUword noflibraries = gt_array_size(r2t->collection);
  GtUword i;
  state->inputfiles_totallength = 0;
  gt_array_sort(r2t->collection, gt_reads2twobit_rli_cmp);
  for (i = 0; i < noflibraries; i++)
  {
    const GtReadsLibraryInfo *rli = gt_array_get(r2t->collection, i);
    state->inputfiles_totallength += rli->total_filelength;
  }
  r2t->twobitencoding = gt_malloc(sizeof (*r2t->twobitencoding) *
    GT_DIVBYUNITSIN2BITENC(state->inputfiles_totallength) + 2UL);
  GT_READS2TWOBIT_INIT_ENCODE_INFO(state->current, r2t->twobitencoding);
  GT_READS2TWOBIT_COPY_ENCODE_INFO(state->current, state->backup);
  for (i = 0; i < (GtUword)(UCHAR_MAX) + 1UL; i++)
    state->char2code[i] = GT_READS2TWOBIT_CODE_UNDEF;
  state->char2code[(unsigned char)'A'] = (GtTwobitencoding)0;
  state->char2code[(unsigned char)'a'] = (GtTwobitencoding)0;
  state->char2code[(unsigned char)'C'] = (GtTwobitencoding)1;
  state->char2code[(unsigned char)'c'] = (GtTwobitencoding)1;
  state->char2code[(unsigned char)'G'] = (GtTwobitencoding)2;
  state->char2code[(unsigned char)'g'] = (GtTwobitencoding)2;
  state->char2code[(unsigned char)'T'] = (GtTwobitencoding)3;
  state->char2code[(unsigned char)'t'] = (GtTwobitencoding)3;
  state->seqlen = 0;
  state->seqlen_mate = 0;
  state->exp_qlen = 0;
  state->varlen_mode = false;
  state->invalid_mode = false;
  state->invalid_sequences = 0;
  state->invalid_total_length = 0;
  gt_assert(r2t->seppos == NULL);
  state->seppos = NULL;
  state->qbuf = NULL;
  state->qbuf_size = 0;
  state->qbuf2 = NULL;
  state->qbuf2_size = 0;
  state->phredbase = r2t->phredbase;
  state->lowqual = r2t->lowqual;
  state->maxlow = r2t->maxlow;
  state->prevcode = GT_READS2TWOBIT_CODE_UNDEF;
  state->hplength = 0;
  state->use_rle = r2t->use_rle;
  state->hplengths = r2t->use_rle ? gt_hplstore_new(
      state->inputfiles_totallength + 2UL) : NULL;
  state->nofh = 0;
  state->hsum = 0;
  state->descs = r2t->descs;
  state->descsfp = r2t->descsfp;
  state->clipdes = r2t->clipdes;
  state->longestdesc = &(r2t->longestdesc);
  state->n_descs = &(r2t->n_descs);
  state->dbuf = gt_str_new();
  state->dbuf2 = gt_str_new();
}

#define GT_READS2TWOBIT_READBUFFER_SIZE ((size_t)256)

#define GT_READS2TWOBIT_WRITECODE_NOCOUNT(EI, CODE, LEN) \
  (LEN)++;\
  (EI).globalpos++;\
  (EI).kmercode = ((EI).kmercode << 2) | (CODE);\
  if (++(EI).codepos == (unsigned short)GT_UNITSIN2BITENC)\
  {\
    *((EI).tbe_next)++ = (EI).kmercode;\
    (EI).codepos = 0;\
    (EI).kmercode = 0;\
  }

#define GT_READS2TWOBIT_WRITECODE(EI, CODE, LEN)\
  GT_READS2TWOBIT_WRITECODE_NOCOUNT(EI, CODE, LEN);\
  (EI).chardistri[CODE]++

static void gt_reads2twobit_init_seppos(GtReads2TwobitEncodeState *state,
    GtUword currentpos)
{
  gt_assert(state->seppos == NULL);
  gt_assert(currentpos > 0);
  gt_assert(state->inputfiles_totallength > currentpos);
  state->seppos_alloc = state->inputfiles_totallength / currentpos;
  gt_log_log("rough estimate of nofseqs = "GT_WU"", state->seppos_alloc);
  state->seppos = gt_malloc(sizeof (state->seppos) * state->seppos_alloc);
  state->current.seppos_nextfree = 0;
}

#define GT_READS2TWOBIT_SEPPOS_INC ((size_t)(1 << 14))

static void gt_reads2twobit_append_seppos(GtReads2TwobitEncodeState *state,
    GtUword pos)
{
  gt_assert(state->seppos != NULL);
  if (state->current.seppos_nextfree == state->seppos_alloc)
  {
    (state->seppos_alloc) += GT_READS2TWOBIT_SEPPOS_INC;
    state->seppos = gt_realloc(state->seppos,
        sizeof (state->seppos) * state->seppos_alloc);
  }
  state->seppos[state->current.seppos_nextfree] = pos;
  state->current.seppos_nextfree++;
}

static void gt_reads2twobit_switch_to_varlen_mode(
    GtReads2TwobitEncodeState *state)
{
  GtUword seqnum;
  GtUword next_seppos;
  gt_assert(state->varlen_mode == false);
  state->varlen_mode = true;
  gt_assert(state->current.nofseqs > 1UL);
  gt_assert(state->current.seqlen_first != state->seqlen);
  gt_assert(state->seqlen > 1UL);
  gt_log_log("readset is varlen: sequences 0.."GT_WU" have length "GT_WU
             ", sequence "GT_WU" is "GT_WU" bp long",
      state->current.nofseqs - 2UL, state->current.seqlen_first - 1UL,
      state->current.nofseqs - 1UL, state->seqlen - 1UL);
  gt_reads2twobit_init_seppos(state,
      state->current.seqlen_first * (state->current.nofseqs - 2UL) +
      state->seqlen);
  next_seppos = 0;
  for (seqnum = 0; seqnum < state->current.nofseqs - 1UL; seqnum++)
  {
    next_seppos += state->current.seqlen_first;
    gt_reads2twobit_append_seppos(state, next_seppos - 1UL);
  }
  gt_assert(next_seppos + state->seqlen == state->current.globalpos);
  gt_reads2twobit_append_seppos(state, state->current.globalpos - 1UL);
  gt_assert(state->current.seppos_nextfree == state->current.nofseqs);
  state->current.seqlen_max = MAX(state->current.seqlen_first, state->seqlen);
  state->current.seqlen_min = MIN(state->current.seqlen_first, state->seqlen);
  state->current.seqlen_first = 0;
}

#define GT_READS2TWOBIT_DEFAULT_SEPARATOR (GtTwobitencoding)3

static inline void gt_reads2twobit_process_sequence_end(
    GtReads2TwobitEncodeState *state)
{
  GT_READS2TWOBIT_WRITECODE_NOCOUNT(state->current,
      GT_READS2TWOBIT_DEFAULT_SEPARATOR, state->seqlen);
  if (state->use_rle)
  {
    gt_assert(state->current.globalpos > 1UL);
    gt_hplstore_set(state->hplengths, state->current.globalpos - 2UL,
        state->hplength);
    state->hsum += state->hplength;
    state->nofh++;
    gt_hplstore_set(state->hplengths, state->current.globalpos - 1UL, 0);
    state->prevcode = GT_READS2TWOBIT_CODE_UNDEF;
  }
  state->exp_qlen++;
  gt_assert(!state->invalid_mode);
  if (state->varlen_mode)
  {
    if (state->seqlen > state->current.seqlen_max)
      state->current.seqlen_max = state->seqlen;
    if (state->seqlen < state->current.seqlen_min)
      state->current.seqlen_min = state->seqlen;
    gt_reads2twobit_append_seppos(state, state->current.globalpos - 1UL);
  }
  else
  {
    if (state->current.nofseqs > 1UL)
    {
      if (state->seqlen != state->current.seqlen_first)
        gt_reads2twobit_switch_to_varlen_mode(state);
    }
    else
    {
      state->current.seqlen_first = state->seqlen;
    }
  }
}

static inline void gt_reads2twobit_prepare_for_new_sequence(
    GtReads2TwobitEncodeState *state)
{
  GT_READS2TWOBIT_COPY_ENCODE_INFO(state->current, state->backup);
  state->current.nofseqs++;
  state->seqlen = 0;
  state->seqlen_mate = 0;
  state->invalid_mode = false;
  state->exp_qlen = 0;
  state->prevcode = GT_READS2TWOBIT_CODE_UNDEF;
  state->hplength = 0;
}

static inline void gt_reads2twobit_process_desc_line(
    GtReads2TwobitEncodeState *state, char *line, bool second_in_pair,
    bool long_line)
{
  GtStr *dbuf;
  if (state->descs == NULL && state->descsfp == NULL)
    return;
  dbuf = second_in_pair ? state->dbuf2 : state->dbuf;
  gt_str_append_cstr(dbuf, long_line ? line : line + 1UL);
}

#define GT_READS2TWOBIT_HANDLE_LONG_DESCRIPTIONS(LINE, FILEPTR, FGETSRETVAL,\
                                                 PROCESS_DESCS, STATE,\
                                                 SECOND_IN_PAIR)\
  /* handle the case in which a description is longer than the line buffer: */\
  while (strlen(LINE) == GT_READS2TWOBIT_READBUFFER_SIZE - (size_t) 1 \
      && (LINE)[GT_READS2TWOBIT_READBUFFER_SIZE - 2] != '\n' \
      && (FGETSRETVAL) != NULL) \
  {\
    (FGETSRETVAL) = \
      fgets((LINE), (int)GT_READS2TWOBIT_READBUFFER_SIZE, (FILEPTR));\
    if ((PROCESS_DESCS) && (FGETSRETVAL) != NULL)\
    {\
      gt_reads2twobit_process_desc_line((STATE), (LINE),\
                                        (SECOND_IN_PAIR), true);\
    }\
  }

static void gt_reads2twobit_switch_to_invalid_mode(
    GtReads2TwobitEncodeState *state)
{
  state->invalid_mode = true;
  state->invalid_sequences++;
  state->invalid_total_length += state->seqlen;
  state->invalid_total_length += state->seqlen_mate;
  if (state->varlen_mode && state->backup.seppos_nextfree == 0)
  {
    gt_log_log("switch back to eqlen mode, seqlen = "GT_WU"",
        state->backup.seqlen_first);
    gt_free(state->seppos);
    state->seppos = NULL;
    state->seppos_alloc = 0;
    state->varlen_mode = false;
  }
  GT_READS2TWOBIT_COPY_ENCODE_INFO(state->backup, state->current);
}

static inline int gt_reads2twobit_process_qualities_line(
    GtReads2TwobitEncodeState *state, const char *line,
    char *qbuf, GtUword *qbuf_next)
{
  GtUword j = 0;
  char c;
  while ((c = line[j++]) != '\0')
  {
    if (c >= state->phredbase)
    {
      if (*qbuf_next == state->exp_qlen)
        return -1;
      qbuf[*qbuf_next] = c - state->phredbase;
      (*qbuf_next)++;
    }
  }
  return 0;
}

static void gt_reads2twobit_apply_quality_filter(
    GtReads2TwobitEncodeState *state, const char *qbuf, GtUword nofq)
{
  GtUword i, low = 0;
  for (i = 0; i < nofq; i++)
    if (qbuf[i] <= state->lowqual)
      low++;
  if (low > state->maxlow)
    gt_reads2twobit_switch_to_invalid_mode(state);
}

static inline void gt_reads2twobit_process_sequence_line(
    GtReads2TwobitEncodeState *state, const char *line)
{
  GtUword j = 0;
  GtTwobitencoding nextcode;
  char c;
  while (true)
  {
    c = line[j++];
    if (!state->invalid_mode && (nextcode = state->char2code[(unsigned char)c])
        != GT_READS2TWOBIT_CODE_UNDEF)
    {
      if (!state->use_rle)
      {
        GT_READS2TWOBIT_WRITECODE(state->current, nextcode, state->seqlen);
      }
      else
      {
        if (nextcode != state->prevcode)
        {
          if (state->seqlen > 0)
          {
            gt_assert(state->current.globalpos > 0);
            gt_hplstore_set(state->hplengths, state->current.globalpos - 1UL,
                state->hplength);
            state->hsum += state->hplength;
            state->nofh++;
            state->hplength = 0;
          }
          GT_READS2TWOBIT_WRITECODE(state->current, nextcode, state->seqlen);
          state->prevcode = nextcode;
        }
        else
        {
          state->hplength += 1;
        }
      }
      state->exp_qlen++;
    }
    else
    {
      if (c == '\0')
        break;
      if (!isspace(c))
      {
        if (!state->invalid_mode)
          gt_reads2twobit_switch_to_invalid_mode(state);
        state->invalid_total_length++;
        state->seqlen++;
        state->exp_qlen++;
      }
    }
  }
}

static int gt_reads2twobit_close_file(FILE *file, GtStr *filename, GtError *err)
{
  int had_err = 0;
  gt_assert(file != NULL);
  if (ferror(file) != 0)
  {
    gt_error_set(err, "Error by reading file %s: %s", gt_str_get(filename),
        strerror(errno));
    had_err = -1;
  }
  gt_fa_fclose(file);
  return had_err;
}

static void gt_reads2twobit_clip_str(GtStr *str)
{
  GtUword i;
  for (i = 0; i < gt_str_length(str); i++)
  {
    if (isspace(gt_str_get(str)[i]))
    {
      gt_str_get(str)[i] = '\n';
      gt_str_set_length(str, i+1);
      return;
    }
  }
}

static void gt_reads2twobit_finalize_descriptions_diskbased(
    GtReads2TwobitEncodeState *state)
{
  GtUword len;
  gt_assert(!state->invalid_mode);
  gt_assert(state->descsfp != NULL);
  gt_assert(gt_str_length(state->dbuf) > 0);
  if (state->clipdes)
    gt_reads2twobit_clip_str(state->dbuf);
  len = gt_str_length(state->dbuf);
  if (len > 0 && gt_str_get(state->dbuf)[len-1] != '\n')
  {
    gt_str_append_char(state->dbuf, '\n');
    len++;
  }
  if (len > (*state->longestdesc))
    (*state->longestdesc) = len;
  gt_xfputs(gt_str_get(state->dbuf), state->descsfp);
  (*state->n_descs)++;
  if (gt_str_length(state->dbuf2) > 0)
  {
    if (state->clipdes)
      gt_reads2twobit_clip_str(state->dbuf2);
    len = gt_str_length(state->dbuf2);
    if (gt_str_get(state->dbuf2)[len-1] != '\n')
    {
      gt_str_append_char(state->dbuf2, '\n');
      len++;
    }
    if (len > (*state->longestdesc))
      (*state->longestdesc) = len;
    gt_xfputs(gt_str_get(state->dbuf2), state->descsfp);
    (*state->n_descs)++;
  }
  gt_str_reset(state->dbuf);
  gt_str_reset(state->dbuf2);
}

static void gt_reads2twobit_finalize_descriptions_membased(
    GtReads2TwobitEncodeState *state)
{
  GtUword i;
  gt_assert(!state->invalid_mode);
  gt_assert(state->descs != NULL);
  gt_assert(gt_str_length(state->dbuf) > 0);
  for (i = 0; i < gt_str_length(state->dbuf); i++)
  {
    gt_desc_buffer_append_char(state->descs, gt_str_get(state->dbuf)[i]);
  }
  gt_desc_buffer_finish(state->descs);
  (*state->n_descs)++;
  if (gt_str_length(state->dbuf2) > 0)
  {
    for (i = 0; i < gt_str_length(state->dbuf2); i++)
    {
      gt_desc_buffer_append_char(state->descs, gt_str_get(state->dbuf2)[i]);
    }
    gt_desc_buffer_finish(state->descs);
    (*state->n_descs)++;
  }
  gt_str_reset(state->dbuf);
  gt_str_reset(state->dbuf2);
}

static void gt_reads2twobit_finalize_descriptions(
    GtReads2TwobitEncodeState *state)
{
  if (state->descs != NULL || state->descsfp != NULL)
  {
    if (state->invalid_mode)
    {
      gt_str_reset(state->dbuf);
      gt_str_reset(state->dbuf2);
    }
    else if (state->descs != NULL)
    {
      gt_assert(state->descsfp == NULL);
      gt_reads2twobit_finalize_descriptions_membased(state);
    }
    else
    {
      gt_assert(state->descsfp != NULL);
      gt_reads2twobit_finalize_descriptions_diskbased(state);
    }
  }
}

static int gt_reads2twobit_encode_unpaired_fastq_library(
    GtReads2TwobitEncodeState *state, GtReadsLibraryInfo *rli, FILE *file,
    char *line)
{
  int had_err = 0;
  GtUword qbuf_next = 0;
  bool qmode = false;
  char *fgetsretval = NULL;
  state->seqlen = 0;
  state->exp_qlen = 0;
  do {
    if (!qmode)
    {
      if (line[0] == '@')
      {
        if (state->current.nofseqs > rli->first_seqnum)
        {
          if (!state->invalid_mode)
            gt_reads2twobit_process_sequence_end(state);
          gt_reads2twobit_finalize_descriptions(state);
        }
        gt_reads2twobit_process_desc_line(state, line, false, false);
        GT_READS2TWOBIT_HANDLE_LONG_DESCRIPTIONS(line, file, fgetsretval,
            true, state, false);
        gt_reads2twobit_prepare_for_new_sequence(state);
      }
      else if (line[0] == '+')
      {
        GT_READS2TWOBIT_HANDLE_LONG_DESCRIPTIONS(line, file, fgetsretval,
            false, state, false);
        if (state->exp_qlen + 1UL > (GtUword)(state->qbuf_size))
        {
          state->qbuf_size = (size_t)state->exp_qlen + 1;
          state->qbuf = gt_realloc(state->qbuf, state->qbuf_size);
        }
        qmode = true;
      }
      else
      {
        gt_reads2twobit_process_sequence_line(state, line);
      }
    }
    else
    {
      had_err = gt_reads2twobit_process_qualities_line(state, line,
          state->qbuf, &qbuf_next);
      if (qbuf_next == state->exp_qlen)
      {
        if (state->maxlow != GT_UNDEF_UWORD)
          gt_reads2twobit_apply_quality_filter(state, state->qbuf, qbuf_next);
        qbuf_next = 0;
        qmode = false;
      }
    }
  } while ((fgetsretval = fgets(line, (int)GT_READS2TWOBIT_READBUFFER_SIZE,
            file)) == line && !had_err);
  gt_reads2twobit_finalize_descriptions(state);
  return had_err;
}

static void gt_reads2twobit_encode_unpaired_fasta_library(
    GtReads2TwobitEncodeState *state, GtReadsLibraryInfo *rli, FILE *file,
    char *line)
{
  char *fgetsretval = NULL;
  do {
    if (line[0] == '>')
    {
      if (state->current.nofseqs > rli->first_seqnum)
      {
        if (!state->invalid_mode)
          gt_reads2twobit_process_sequence_end(state);
        gt_reads2twobit_finalize_descriptions(state);
      }
      gt_reads2twobit_process_desc_line(state, line, false, false);
      GT_READS2TWOBIT_HANDLE_LONG_DESCRIPTIONS(line, file, fgetsretval,
          true, state, false);
      gt_reads2twobit_prepare_for_new_sequence(state);
    }
    else if (!state->invalid_mode)
    {
      gt_reads2twobit_process_sequence_line(state, line);
    }
  } while ((fgetsretval = fgets(line, (int)GT_READS2TWOBIT_READBUFFER_SIZE,
            file)) == line);
  gt_reads2twobit_finalize_descriptions(state);
}

static int gt_reads2twobit_encode_unpaired_library(
    GtReads2TwobitEncodeState *state, GtReadsLibraryInfo *rli, GtError *err)
{
  int had_err = 0;
  FILE *file;
  file = gt_fa_fopen(gt_str_get(rli->filename1), "r", err);
  if (file == NULL)
    had_err = -1;
  if (!had_err)
  {
    const GtUword invalid_tl_before = state->invalid_total_length,
          invalid_s_before = state->invalid_sequences;
    char line[GT_READS2TWOBIT_READBUFFER_SIZE], *fgetsretval;
    rli->first_seqnum = state->current.nofseqs;
    fgetsretval = fgets(line, (int)GT_READS2TWOBIT_READBUFFER_SIZE, file);
    if (fgetsretval == line)
    {
      if (line[0] == '>')
      {
        gt_reads2twobit_encode_unpaired_fasta_library(state, rli, file, line);
      }
      else if (line[0] == '@')
      {
        had_err = gt_reads2twobit_encode_unpaired_fastq_library(state, rli,
            file, line);
        if (had_err != 0)
          gt_error_set(err, "%s: error in FASTQ format",
              gt_str_get(rli->filename1));
      }
      else
      {
        gt_error_set(err, "%s: unknown format", gt_str_get(rli->filename1));
        had_err = -1;
      }
    }
    else
    {
      /* empty file or error */
      gt_assert(fgetsretval == NULL);
    }
    if (!state->invalid_mode)
      gt_reads2twobit_process_sequence_end(state);
    gt_assert(state->current.nofseqs >= rli->first_seqnum);
    rli->nofseqs = state->current.nofseqs - rli->first_seqnum;
    rli->total_seqlength = state->varlen_mode
      ?  state->seppos[state->current.nofseqs - 1UL] + 1UL -
         ((rli->first_seqnum == 0) ? 0 : state->seppos[rli->first_seqnum - 1UL])
      : state->current.seqlen_first * rli->nofseqs;
    /* the following is not necessary, but is useful for the tests */
    rli->total_filelength -= (state->invalid_total_length - invalid_tl_before +
        3UL * (state->invalid_sequences - invalid_s_before));
    if (!had_err)
      had_err = gt_reads2twobit_close_file(file, rli->filename1, err);
    else
      gt_fa_fclose(file);
  }
  return had_err;
}

static inline int gt_reads2twobit_process_fastq_mate_pair(
    GtReads2TwobitEncodeState *state, char *line2, FILE *file2, bool *file2new)
{
  int had_err = 0;
  char *fgetsretval = NULL;
  bool was_invalid = state->invalid_mode;
  GtUword prev_seqlen = state->seqlen;
  GtUword qbuf2_next = 0;
  bool qmode = false;
  gt_assert((state->descs == NULL && state->descsfp == NULL) ||
             gt_str_length(state->dbuf) > 0);
  if (*file2new)
  {
    fgetsretval = fgets(line2, (int)GT_READS2TWOBIT_READBUFFER_SIZE, file2);
    gt_assert(fgetsretval == line2);
    *file2new = false;
  }
  /* at this point one line2 is always the (beginnning of the) @ description */
  gt_assert(line2[0] == '@');
  gt_reads2twobit_process_desc_line(state, line2, true, false);
  GT_READS2TWOBIT_HANDLE_LONG_DESCRIPTIONS(line2, file2, fgetsretval, true,
      state, true);
  gt_assert((state->descs == NULL && state->descsfp == NULL) ||
             gt_str_length(state->dbuf2) > 0);
  state->seqlen_mate = state->seqlen;
  state->seqlen = 0;
  state->exp_qlen = 0;
  if (!state->invalid_mode)
    state->current.nofseqs++;
  else
    state->invalid_sequences++;
  while (!had_err)
  {
    fgetsretval = fgets(line2, (int)GT_READS2TWOBIT_READBUFFER_SIZE, file2);
    if (fgetsretval == NULL)
      break;
    if (!qmode)
    {
      if (line2[0] == '@')
        break;
      else if (line2[0] == '+')
      {
        GT_READS2TWOBIT_HANDLE_LONG_DESCRIPTIONS(line2, file2, fgetsretval,
            false, state, true);
        if (state->exp_qlen + 1UL > (GtUword)(state->qbuf2_size))
        {
          state->qbuf2_size = (size_t)state->exp_qlen + 1;
          state->qbuf2 = gt_realloc(state->qbuf2, state->qbuf2_size);
        }
        qmode = true;
      }
      else
      {
        gt_reads2twobit_process_sequence_line(state, line2);
      }
    }
    else
    {
      had_err = gt_reads2twobit_process_qualities_line(state, line2,
          state->qbuf2, &qbuf2_next);
      if (had_err == -1)
        had_err = -2;
      if (!had_err && qbuf2_next == state->exp_qlen)
      {
        if (state->maxlow != GT_UNDEF_UWORD)
          gt_reads2twobit_apply_quality_filter(state, state->qbuf2, qbuf2_next);
        qbuf2_next = 0;
        qmode = false;
      }
    }
  }
  if (!had_err)
    gt_reads2twobit_finalize_descriptions(state);
  if (!had_err && !state->invalid_mode)
    gt_reads2twobit_process_sequence_end(state);
  if (!had_err && !was_invalid && state->invalid_mode)
  {
    state->invalid_sequences++;
    state->invalid_total_length += (prev_seqlen - 1UL);
  }
  return had_err;
}

static int gt_reads2twobit_encode_interleaved_paired_fastq_library(
    GtReads2TwobitEncodeState *state, GtReadsLibraryInfo *rli, FILE *file,
    char *line)
{
  int had_err = 0;
  GtUword qbuf_next = 0;
  bool qmode = false;
  bool processing_mate = true;
  char *fgetsretval = line;
  state->seqlen = 0;
  state->exp_qlen = 0;
  do {
    if (!qmode)
    {
      if (line[0] == '@')
      {
        processing_mate = !processing_mate;
        if (state->current.nofseqs > rli->first_seqnum)
        {
          if (!state->invalid_mode)
            gt_reads2twobit_process_sequence_end(state);
        }
        if (!processing_mate && (state->current.nofseqs > rli->first_seqnum ||
              state->invalid_mode))
        {
          gt_reads2twobit_finalize_descriptions(state);
        }
        gt_reads2twobit_process_desc_line(state, line, processing_mate,
            false);
        GT_READS2TWOBIT_HANDLE_LONG_DESCRIPTIONS(line, file, fgetsretval,
            true, state, processing_mate);
        if (processing_mate)
        {
          if (!state->invalid_mode)
            state->current.nofseqs++;
          state->seqlen_mate = state->seqlen;
          state->seqlen = 0;
          state->exp_qlen = 0;
        }
        else
        {
          gt_reads2twobit_prepare_for_new_sequence(state);
        }
      }
      else if (line[0] == '+')
      {
        GT_READS2TWOBIT_HANDLE_LONG_DESCRIPTIONS(line, file, fgetsretval,
            false, state, processing_mate);
        if (state->exp_qlen + 1UL > (GtUword)(state->qbuf_size))
        {
          state->qbuf_size = (size_t)state->exp_qlen + 1;
          state->qbuf = gt_realloc(state->qbuf, state->qbuf_size);
        }
        qmode = true;
      }
      else
      {
        gt_reads2twobit_process_sequence_line(state, line);
      }
    }
    else
    {
      had_err = gt_reads2twobit_process_qualities_line(state, line,
          state->qbuf, &qbuf_next);
      if (qbuf_next == state->exp_qlen)
      {
        if (state->maxlow != GT_UNDEF_UWORD)
          gt_reads2twobit_apply_quality_filter(state, state->qbuf, qbuf_next);
        qbuf_next = 0;
        qmode = false;
      }
    }
  } while ((fgetsretval = fgets(line, (int)GT_READS2TWOBIT_READBUFFER_SIZE,
            file)) == line && !had_err);
  if (!had_err && !state->invalid_mode)
    gt_reads2twobit_process_sequence_end(state);
  if (state->current.nofseqs > rli->first_seqnum)
    gt_reads2twobit_finalize_descriptions(state);
  return had_err;
}

static int gt_reads2twobit_encode_twofile_paired_fastq_library(
    GtReads2TwobitEncodeState *state, GtReadsLibraryInfo *rli, FILE *file1,
    FILE *file2, char *line1)
{
  int had_err = 0;
  GtUword qbuf_next = 0;
  bool qmode = false;
  bool file2new = true;
  char *fgetsretval = NULL,
       line2[GT_READS2TWOBIT_READBUFFER_SIZE];
  line2[0] = '\0';
  state->seqlen = 0;
  state->exp_qlen = 0;
  do {
    if (!qmode)
    {
      if (line1[0] == '@')
      {
        if (state->current.nofseqs > rli->first_seqnum)
        {
          if (!state->invalid_mode)
            gt_reads2twobit_process_sequence_end(state);
          had_err = gt_reads2twobit_process_fastq_mate_pair(state, line2,
              file2, &file2new);
        }
        else
        {
          if (state->invalid_mode)
            had_err = gt_reads2twobit_process_fastq_mate_pair(state, line2,
                file2, &file2new);
        }
        gt_reads2twobit_process_desc_line(state, line1, false, false);
        GT_READS2TWOBIT_HANDLE_LONG_DESCRIPTIONS(line1, file1, fgetsretval,
            true, state, false);
        gt_assert((state->descs == NULL && state->descsfp == NULL) ||
                   gt_str_length(state->dbuf) > 0);
        gt_reads2twobit_prepare_for_new_sequence(state);
      }
      else if (line1[0] == '+')
      {
        GT_READS2TWOBIT_HANDLE_LONG_DESCRIPTIONS(line1, file1, fgetsretval,
            false, state, false);
        if (state->exp_qlen + 1UL > (GtUword)(state->qbuf_size))
        {
          state->qbuf_size = (size_t)state->exp_qlen + 1;
          state->qbuf = gt_realloc(state->qbuf, state->qbuf_size);
        }
        qmode = true;
      }
      else
      {
        gt_reads2twobit_process_sequence_line(state, line1);
      }
    }
    else
    {
      had_err = gt_reads2twobit_process_qualities_line(state, line1,
          state->qbuf, &qbuf_next);
      if (qbuf_next == state->exp_qlen)
      {
        if (state->maxlow != GT_UNDEF_UWORD)
          gt_reads2twobit_apply_quality_filter(state, state->qbuf, qbuf_next);
        qbuf_next = 0;
        qmode = false;
      }
    }
  } while ((fgetsretval = fgets(line1, (int)GT_READS2TWOBIT_READBUFFER_SIZE,
            file1)) == line1 && !had_err);
  if (!had_err && !state->invalid_mode)
    gt_reads2twobit_process_sequence_end(state);
  if (!had_err)
  {
    had_err = gt_reads2twobit_process_fastq_mate_pair(state, line2, file2,
        &file2new);
  }
  return had_err;
}

static inline void gt_reads2twobit_process_fasta_mate_pair(
    GtReads2TwobitEncodeState *state, char *line2, FILE *file2)
{
  char *fgetsretval = NULL;
  bool was_invalid = state->invalid_mode;
  GtUword prev_seqlen = state->seqlen;
  if (line2[0] != '>')
  {
    fgetsretval = fgets(line2, (int)GT_READS2TWOBIT_READBUFFER_SIZE,
        file2);
    gt_assert(fgetsretval == line2);
  }
  gt_assert(line2[0] == '>');
  gt_reads2twobit_process_desc_line(state, line2, true, false);
  GT_READS2TWOBIT_HANDLE_LONG_DESCRIPTIONS(line2, file2, fgetsretval, true,
      state, true);
  state->seqlen = 0;
  state->exp_qlen = 0;
  state->seqlen_mate = state->seqlen;
  if (!state->invalid_mode)
    state->current.nofseqs++;
  else
    state->invalid_sequences++;
  while (true)
  {
    fgetsretval = fgets(line2, (int)GT_READS2TWOBIT_READBUFFER_SIZE, file2);
    if (fgetsretval != line2 || line2[0] == '>')
      break;
    gt_reads2twobit_process_sequence_line(state, line2);
  }
  gt_reads2twobit_finalize_descriptions(state);
  if (!state->invalid_mode)
    gt_reads2twobit_process_sequence_end(state);
  if (!was_invalid && state->invalid_mode)
  {
    state->invalid_sequences++;
    state->invalid_total_length += (prev_seqlen - 1UL);
  }
}

static void gt_reads2twobit_encode_interleaved_paired_fasta_library(
    GtReads2TwobitEncodeState *state, GtReadsLibraryInfo *rli, FILE *file,
    char *line)
{
  char *fgetsretval = NULL;
  bool processing_mate = true;
  do {
    if (line[0] == '>')
    {
      processing_mate = !processing_mate;
      if (state->current.nofseqs > rli->first_seqnum)
      {
        if (!state->invalid_mode)
          gt_reads2twobit_process_sequence_end(state);
        if (!processing_mate)
          gt_reads2twobit_finalize_descriptions(state);
      }
      gt_reads2twobit_process_desc_line(state, line, processing_mate, false);
      GT_READS2TWOBIT_HANDLE_LONG_DESCRIPTIONS(line, file, fgetsretval, true,
          state, processing_mate);
      if (processing_mate)
      {
        if (!state->invalid_mode)
          state->current.nofseqs++;
        state->seqlen_mate = state->seqlen;
        state->seqlen = 0;
        state->exp_qlen = 0;
      }
      else
      {
        gt_reads2twobit_prepare_for_new_sequence(state);
      }
    }
    else
    {
      gt_reads2twobit_process_sequence_line(state, line);
    }
  } while ((fgetsretval = fgets(line, (int)GT_READS2TWOBIT_READBUFFER_SIZE,
            file)) == line);
  if (!state->invalid_mode)
    gt_reads2twobit_process_sequence_end(state);
  gt_reads2twobit_finalize_descriptions(state);
}

static void gt_reads2twobit_encode_twofile_paired_fasta_library(
    GtReads2TwobitEncodeState *state, GtReadsLibraryInfo *rli, FILE *file1,
    FILE *file2, char *line1)
{
  char *fgetsretval = NULL,
       line2[GT_READS2TWOBIT_READBUFFER_SIZE];
  line2[0] = '\0';
  do {
    if (line1[0] == '>')
    {
      if (state->current.nofseqs > rli->first_seqnum)
      {
        if (!state->invalid_mode)
          gt_reads2twobit_process_sequence_end(state);
        gt_reads2twobit_process_fasta_mate_pair(state, line2, file2);
      }
      else
      {
        if (state->invalid_mode)
          gt_reads2twobit_process_fasta_mate_pair(state, line2, file2);
      }
      gt_reads2twobit_process_desc_line(state, line1, false, false);
      GT_READS2TWOBIT_HANDLE_LONG_DESCRIPTIONS(line1, file1, fgetsretval,
          true, state, false);
      gt_reads2twobit_prepare_for_new_sequence(state);
    }
    else if (!state->invalid_mode)
    {
      gt_reads2twobit_process_sequence_line(state, line1);
    }
  } while ((fgetsretval = fgets(line1, (int)GT_READS2TWOBIT_READBUFFER_SIZE,
            file1)) == line1);
  if (!state->invalid_mode)
    gt_reads2twobit_process_sequence_end(state);
  gt_reads2twobit_process_fasta_mate_pair(state, line2, file2);
}

static int gt_reads2twobit_encode_paired_library(
    GtReads2TwobitEncodeState *state, GtReadsLibraryInfo *rli, GtError *err)
{
  int had_err = 0;
  FILE *file1, *file2 = NULL;
  gt_assert(rli->filename1 != NULL);
  file1 = gt_fa_fopen(gt_str_get(rli->filename1), "r", err);
  if (file1 == NULL)
    had_err = -1;
  if (!had_err && rli->filename2 != NULL)
  {
    file2 = gt_fa_fopen(gt_str_get(rli->filename2), "r", err);
    if (file2 == NULL)
      had_err = -1;
  }
  if (!had_err)
  {
    const GtUword invalid_tl_before = state->invalid_total_length,
          invalid_s_before = state->invalid_sequences;
    char line1[GT_READS2TWOBIT_READBUFFER_SIZE],
         *fgetsretval;
    rli->first_seqnum = state->current.nofseqs;
    fgetsretval = fgets(line1, (int)GT_READS2TWOBIT_READBUFFER_SIZE, file1);
    if (fgetsretval == line1)
    {
      if (line1[0] == '>')
      {
        if (file2 == NULL)
        {
          gt_log_log("encode interleaved paired fasta library");
          gt_reads2twobit_encode_interleaved_paired_fasta_library(state, rli,
              file1, line1);
        }
        else
        {
          gt_log_log("encode two-file paired fasta library");
          gt_reads2twobit_encode_twofile_paired_fasta_library(state, rli, file1,
              file2, line1);
        }
      }
      else if (line1[0] == '@')
      {
        if (file2 == NULL)
        {
          gt_log_log("encode interleaved paired fastQ library");
          had_err = gt_reads2twobit_encode_interleaved_paired_fastq_library(
              state, rli, file1, line1);
        }
        else
        {
          gt_log_log("encode two-file paired fastQ library");
          had_err = gt_reads2twobit_encode_twofile_paired_fastq_library(
              state, rli, file1, file2, line1);
        }
        if (had_err != 0)
          gt_error_set(err, "%s: error in FASTQ format", had_err == -1 ?
              gt_str_get(rli->filename1) : gt_str_get(rli->filename2));
      }
      else
      {
        gt_error_set(err, "%s: unknown format", gt_str_get(rli->filename1));
        had_err = -1;
      }
    }
    else
    {
      /* empty file or error */
      gt_assert(fgetsretval == NULL);
    }
    gt_assert(state->current.nofseqs >= rli->first_seqnum);
    rli->nofseqs = state->current.nofseqs - rli->first_seqnum;
    rli->total_seqlength = state->varlen_mode
      ?  state->seppos[state->current.nofseqs - 1UL] + 1UL -
         ((rli->first_seqnum == 0) ? 0 : state->seppos[rli->first_seqnum - 1UL])
      : state->current.seqlen_first * rli->nofseqs;
    /* the following is not necessary, but is useful for the tests */
    rli->total_filelength -= (state->invalid_total_length - invalid_tl_before +
        3UL * (state->invalid_sequences - invalid_s_before));
  }
  if (!had_err)
    had_err = gt_reads2twobit_close_file(file1, rli->filename1, err);
  else if (file1 != NULL)
    gt_fa_fclose(file1);
  if (!had_err && file2 != NULL)
    had_err = gt_reads2twobit_close_file(file2, rli->filename2, err);
  else if (file2 != NULL)
    gt_fa_fclose(file2);
  return had_err;
}

static void gt_reads2twobit_tbe_flush_and_realloc(GtReads2Twobit *r2t,
    GtReads2TwobitEncodeState *state)
{
  if (state->current.codepos > 0)
  {
    GtUword shift = (GT_UNITSIN2BITENC - state->current.codepos) << 1UL;
    *(state->current.tbe_next++) =
      (GtTwobitencoding)(state->current.kmercode << shift);
  }
  if (state->current.nofseqs > 0)
  {
    gt_log_log("realloc tbe, total_seqlength = "GT_WU", nofseqs = "GT_WU"",
        r2t->total_seqlength, state->current.nofseqs);
    gt_assert(r2t->total_seqlength > 0);
    r2t->twobitencoding = gt_realloc(r2t->twobitencoding,
        sizeof (*r2t->twobitencoding) *
        (GT_DIVBYUNITSIN2BITENC(r2t->total_seqlength) + 2UL));
  }
  else
  {
    gt_free(r2t->twobitencoding);
    r2t->twobitencoding = NULL;
  }
}

static void gt_reads2twobit_finalize_encode(GtReads2Twobit *r2t,
    GtReads2TwobitEncodeState *state)
{
  unsigned int i;
  r2t->nofseqs = state->current.nofseqs;
  r2t->current_sepcode = GT_READS2TWOBIT_DEFAULT_SEPARATOR;
  r2t->invalid_sequences = state->invalid_sequences;
  r2t->invalid_total_length = state->invalid_total_length;
  for (i = 0; i < GT_READS2TWOBIT_ALPHASIZE; i++)
    r2t->chardistri[i] = state->current.chardistri[i];
  if (state->varlen_mode)
  {
    state->seppos = gt_realloc(state->seppos, sizeof (state->seppos) *
        (state->current.nofseqs));
    r2t->seqlen_eqlen = 0;
    r2t->seqlen_max = state->current.seqlen_max;
    r2t->seqlen_min = state->current.seqlen_min;
    r2t->total_seqlength = state->seppos[state->current.nofseqs - 1UL];
    r2t->seppos = state->seppos;
  }
  else
  {
    r2t->seqlen_eqlen = state->current.seqlen_first;
    r2t->seqlen_max = state->current.seqlen_first;
    r2t->seqlen_min = state->current.seqlen_first;
    if (state->current.seqlen_first > 0)
      r2t->total_seqlength = state->current.seqlen_first *
        state->current.nofseqs - 1UL;
    else
      r2t->total_seqlength = 0;
    gt_assert(state->seppos == NULL);
  }
  if (r2t->use_rle)
  {
    if (r2t->total_seqlength > 0)
    {
      r2t->hplengths = state->hplengths;
      gt_hplstore_finalize(r2t->hplengths, r2t->total_seqlength);
      r2t->approx_avhlen = (double)state->hsum / state->nofh + 1;
    }
    else
    {
      gt_hplstore_delete(state->hplengths);
    }
  }
  if (state->qbuf_size > 0)
    gt_free(state->qbuf);
  if (state->qbuf2_size > 0)
    gt_free(state->qbuf2);
  gt_str_delete(state->dbuf);
  gt_str_delete(state->dbuf2);
  gt_reads2twobit_tbe_flush_and_realloc(r2t, state);
}

int gt_reads2twobit_encode(GtReads2Twobit *r2t, GtError *err)
{
  int had_err = 0;
  const GtUword noflibraries = gt_array_size(r2t->collection);
  GtUword libnum;
  GtReads2TwobitEncodeState state;

  gt_error_check(err);
  gt_assert(r2t != NULL);
  gt_assert(r2t->twobitencoding == NULL);
  gt_reads2twobit_init_encode(r2t, &state);
  for (libnum = 0; libnum < noflibraries && !had_err; libnum++)
  {
    GtReadsLibraryInfo *rli = gt_array_get(r2t->collection, libnum);
    if (!rli->paired)
    {
      had_err = gt_reads2twobit_encode_unpaired_library(&state, rli, err);
    }
    else
    {
      had_err = gt_reads2twobit_encode_paired_library(&state, rli, err);
    }
  }
  gt_reads2twobit_finalize_encode(r2t, &state);
  return had_err;
}

static void gt_reads2twobit_collect_fileinfo(const GtReads2Twobit *r2t,
    GtFilelengthvalues **filelengthtab, GtStrArray **filenametab)
{
  GtUword i, noflibraries;
  noflibraries = gt_array_size(r2t->collection);
  *filenametab = gt_str_array_new();
  if (filelengthtab != NULL)
  {
    *filelengthtab = gt_malloc(sizeof (**filelengthtab) * noflibraries);
  }
  for (i = 0; i < noflibraries; i++)
  {
    GtStr *libname;
    GtReadsLibraryInfo *rli;
    rli = gt_array_get(r2t->collection, i);
    if (filelengthtab != NULL)
    {
      (*filelengthtab)[i].effectivelength =
        (uint64_t)rli->total_seqlength - 1UL;
      if (r2t->seqlen_eqlen == 0 && i == noflibraries - 1UL)
        (*filelengthtab)[i].effectivelength--;
      (*filelengthtab)[i].length = (uint64_t)rli->total_filelength;
    }
    libname = gt_str_clone(rli->filename1);
    if (rli->filename2 != NULL)
    {
      gt_str_append_char(libname, GT_READS2TWOBIT_LIBSPECSEP);
      gt_str_append_str(libname, rli->filename2);
    }
    if (rli->paired)
    {
      gt_str_append_char(libname, GT_READS2TWOBIT_LIBSPECSEP);
      gt_str_append_uword(libname, rli->insertlength);
      if (rli->stdev > 0)
      {
        gt_str_append_char(libname, GT_READS2TWOBIT_INSERTSEP);
        gt_str_append_uword(libname, rli->stdev);
      }
    }
    gt_str_array_add(*filenametab, libname);
    gt_str_delete(libname);
  }
}

static inline GtTwobitencoding gt_reads2twobit_less_frequent_char(
    GtReads2Twobit *r2t)
{
  GtTwobitencoding i, code;
  GtUword lowest_value, value;
  lowest_value = r2t->chardistri[0];
  code = 0;
  for (i = (GtTwobitencoding)1;
      i < (GtTwobitencoding)GT_READS2TWOBIT_ALPHASIZE; i++)
  {
    value = r2t->chardistri[i];
    if (value < lowest_value)
    {
      lowest_value = value;
      code = i;
    }
  }
  gt_log_log("less frequent char code: "GT_WU"", (GtUword)code);
  return code;
}

static void gt_reads2twobit_zeropad_tbe(GtReads2Twobit *r2t)
{
  GtUword pos, codenum, posincode, shift;
  pos = r2t->total_seqlength - 1UL;
  codenum = GT_DIVBYUNITSIN2BITENC(pos);
  posincode = GT_MODBYUNITSIN2BITENC(pos);
  if (posincode < (GtUword)GT_UNITSIN2BITENC - 1UL)
  {
    shift = GT_MULT2(GT_UNITSIN2BITENC - 1UL - posincode);
    r2t->twobitencoding[codenum] =
      (GtTwobitencoding)((r2t->twobitencoding[codenum] >> shift) << shift);
  }
  r2t->twobitencoding[codenum + 1UL] = 0;
}

static void gt_reads2twobit_seek_sequence(const GtReads2Twobit *r2t,
    GtUword seqnum, GtUword *seqlen, GtTwobitencoding *firstcode,
    GtUword *charsinfirstcode, GtTwobitencoding **nextcode_ptr)
{
  GtUword pos;
  *seqlen = r2t->seqlen_eqlen;
  if (*seqlen == 0)
  {
    if (seqnum == 0)
    {
      *seqlen = r2t->seppos[0] + 1UL;
      pos = 0;
    }
    else
    {
      *seqlen = r2t->seppos[seqnum] - r2t->seppos[seqnum - 1UL];
      pos = r2t->seppos[seqnum - 1UL] + 1UL;
    }
  }
  else
  {
    pos = seqnum * (*seqlen);
  }
  *nextcode_ptr = r2t->twobitencoding + GT_DIVBYUNITSIN2BITENC(pos);
  *charsinfirstcode = GT_UNITSIN2BITENC - GT_MODBYUNITSIN2BITENC(pos);
  *firstcode = *((*nextcode_ptr)++);
}

/* decodes the specified sequence in Fasta format; the <decoded> buffer
   must be large enough */
void gt_reads2twobit_decode_sequence(const GtReads2Twobit *r2t,
    GtUword seqnum, char *decoded)
{
  GtTwobitencoding code;
  GtUword pos, seqlen, charsincode;
  GtTwobitencoding *nextencoded;
  const char code2char[] = "acgt";
  char *nextdecoded = decoded;

  gt_reads2twobit_seek_sequence(r2t, seqnum, &seqlen, &code, &charsincode,
      &nextencoded);
  *(nextdecoded++) = '>';
  *(nextdecoded++) = '\n';
  for (pos = 0; pos < seqlen - 1UL; pos++)
  {
    if (charsincode == 0)
    {
      code = *(nextencoded++);
      charsincode = (GtUword)GT_UNITSIN2BITENC;
    }
    *(nextdecoded++) = code2char[code >> ((--charsincode) << 1) & 3];
  }
  *(nextdecoded++) = '\n';
  *(nextdecoded) = '\0';
}

static GtUword gt_reads2twobit_subtract_from_chardistri(
    GtReads2Twobit *r2t, GtUword seqnum)
{
  GtTwobitencoding code;
  GtUword pos, seqlen, charsincode;
  GtTwobitencoding *nextencoded;

  gt_reads2twobit_seek_sequence(r2t, seqnum, &seqlen, &code, &charsincode,
      &nextencoded);
  for (pos = 0; pos < seqlen - 1UL; pos++)
  {
    if (charsincode == 0)
    {
      code = *(nextencoded++);
      charsincode = (GtUword)GT_UNITSIN2BITENC;
    }
    r2t->chardistri[code >> ((--charsincode) << 1) & 3]--;
  }
  return seqlen;
}

/* decodes the sequences <seqnum_from> to <seqnum_from>+<nofseqs>-1
   in MultiFasta format and outputs to <outfp>; if <skip> is not NULL,
   then skips any sequence for which the corresponding bit is set */
void gt_reads2twobit_decode_range(const GtReads2Twobit *r2t,
    GtFile *outfp, GtUword seqnum_from, GtUword nofseqs,
    const GtBitsequence *skip)
{
  GtTwobitencoding code;
  unsigned short charsincode;
  const char code2char[] = "acgt";
  GtUword seqnum, pos, nextsep, nextdecoded, seqnum_to;
  const GtTwobitencoding *nextencoded;
  char *decoded;

  gt_assert(r2t->seqlen_max > 0);
  if (nofseqs == 0)
    return;
  seqnum_to = seqnum_from + nofseqs - 1UL;
  decoded = gt_malloc(sizeof (*decoded) * (r2t->seqlen_max + 3UL));
  decoded[0] = '>';
  decoded[1] = '\n';
  nextdecoded = 2UL;

  seqnum = seqnum_from;
  if (skip != NULL)
    while (GT_ISIBITSET(skip, seqnum))
      seqnum++;

  if (r2t->seqlen_eqlen > 0)
  {
    nextsep = r2t->seqlen_eqlen * (seqnum + 1UL) - 1UL;
    pos = seqnum * r2t->seqlen_eqlen;
  }
  else
  {
    nextsep = r2t->seppos[seqnum];
    pos = (seqnum == 0) ? 0 : r2t->seppos[seqnum - 1UL] + 1UL;
  }

  nextencoded = r2t->twobitencoding + GT_DIVBYUNITSIN2BITENC(pos);
  code = *(nextencoded++);
  charsincode = (unsigned short)GT_UNITSIN2BITENC -
    (unsigned short)GT_MODBYUNITSIN2BITENC(pos);

  while (true)
  {
    if (charsincode == 0)
    {
      code = *(nextencoded++);
      charsincode = (unsigned short)GT_UNITSIN2BITENC;
    }
    if (pos++ < nextsep)
    {
      gt_assert(nextsep - pos <= r2t->seqlen_max);
      decoded[nextdecoded++] = code2char[code >> ((--charsincode) << 1) & 3];
    }
    else
    {
      /* output sequence */
      charsincode--;
      decoded[nextdecoded++] = '\n';
      decoded[nextdecoded] = '\0';
      gt_file_xfputs(decoded, outfp);
      nextdecoded = 2UL;

      seqnum++;
      if (seqnum > seqnum_to)
        break;
      if (skip != NULL && GT_ISIBITSET(skip, seqnum))
      {
        /* jump to next non-contained sequence */
        while (seqnum <= seqnum_to && GT_ISIBITSET(skip, seqnum))
          seqnum++;
        if (seqnum > seqnum_to)
          break;
        pos = r2t->seqlen_eqlen > 0 ? seqnum * r2t->seqlen_eqlen :
          r2t->seppos[seqnum - 1UL] + 1UL;
        nextencoded = r2t->twobitencoding + GT_DIVBYUNITSIN2BITENC(pos);
        code = *(nextencoded++);
        charsincode = (unsigned short)GT_UNITSIN2BITENC -
          (unsigned short)GT_MODBYUNITSIN2BITENC(pos);
      }
      nextsep = (r2t->seqlen_eqlen > 0)
        ? r2t->seqlen_eqlen * (seqnum + 1UL) - 1UL
        : r2t->seppos[seqnum];
    }
  }
  gt_free(decoded);
}

int gt_reads2twobit_write_fasta(const GtReads2Twobit *r2t, char *path,
    GtBitsequence *skip, GtError *err)
{
  GtFile *file;
  file = gt_file_new(path, "w", err);
  if (file == NULL)
    return -1;
  gt_reads2twobit_decode_range(r2t, file, 0, r2t->nofseqs, skip);
  gt_file_delete(file);
  return 0;
}

static void gt_reads2twobit_write_encoded_nocodesshift(GtReads2Twobit *r2t,
    GtTwobitencoding *outputbuffer, GtTwobitencoding offset,
    GtUword firstcodeidx, GtUword lastcodeidx)
{
  if (offset == 0)
    outputbuffer[0] = r2t->twobitencoding[firstcodeidx];
  else
  {
    GtTwobitencoding mask;
    mask =
      ((GtTwobitencoding)1 << GT_MULT2(GT_UNITSIN2BITENC - offset)) - 1;
    outputbuffer[0] =
      (r2t->twobitencoding[firstcodeidx] & mask) |
      (outputbuffer[0] & ~mask);
  }
  if (lastcodeidx > firstcodeidx)
    memcpy(outputbuffer + 1, r2t->twobitencoding + firstcodeidx + 1,
        sizeof (GtTwobitencoding) * (lastcodeidx - firstcodeidx));
}

static void gt_reads2twobit_write_encoded_leftcodesshift(
    GtReads2Twobit *r2t, GtTwobitencoding *outputbuffer,
    GtTwobitencoding inputoffset, GtTwobitencoding outputoffset,
    GtUword firstcodeidx, GtUword lastcodeidx)
{
  const GtTwobitencoding netoffset = inputoffset - outputoffset;
  const GtTwobitencoding shiftright =
    GT_MULT2(GT_UNITSIN2BITENC - netoffset);
  const GtTwobitencoding shiftleft = GT_MULT2(netoffset);
  GtUword i;
  GtTwobitencoding *nextinoutputbuffer = outputbuffer;
  if (outputoffset == 0)
    outputbuffer[0] = r2t->twobitencoding[firstcodeidx] << shiftleft;
  else
  {
    GtTwobitencoding mask;
    mask = ((GtTwobitencoding)1 << GT_MULT2(GT_UNITSIN2BITENC -
          outputoffset)) - 1;
    outputbuffer[0] =
      ((r2t->twobitencoding[firstcodeidx] << shiftleft) & mask) |
      (outputbuffer[0] & ~mask);
  }
  for (i = firstcodeidx + 1UL; i <= lastcodeidx; i++)
  {
    *(nextinoutputbuffer) |= (r2t->twobitencoding[i] >> shiftright);
    *(++nextinoutputbuffer) = (r2t->twobitencoding[i] << shiftleft);
  }
}

static void gt_reads2twobit_write_encoded_rightcodesshift(
    GtReads2Twobit *r2t, GtTwobitencoding *outputbuffer,
    GtTwobitencoding inputoffset, GtTwobitencoding outputoffset,
    GtUword firstcodeidx, GtUword lastcodeidx)
{
  const GtTwobitencoding netoffset = outputoffset - inputoffset;
  const GtTwobitencoding shiftright = GT_MULT2(netoffset);
  const GtTwobitencoding shiftleft =
    GT_MULT2(GT_UNITSIN2BITENC - netoffset);
  GtUword i;
  GtTwobitencoding *nextinoutputbuffer = outputbuffer + 1;
  GtTwobitencoding mask;
  mask = ((GtTwobitencoding)1 << GT_MULT2(GT_UNITSIN2BITENC -
        outputoffset)) - 1;
  outputbuffer[0] =
    ((r2t->twobitencoding[firstcodeidx] >> shiftright) & mask) |
    (outputbuffer[0] & ~mask);
  outputbuffer[1] =
    r2t->twobitencoding[firstcodeidx] << shiftleft;
  for (i = firstcodeidx + 1; i <= lastcodeidx; i++)
  {
    *(nextinoutputbuffer) |= (r2t->twobitencoding[i] >> shiftright);
    *(++nextinoutputbuffer) = (r2t->twobitencoding[i] << shiftleft);
  }
}

GtTwobitencoding* gt_reads2twobit_write_encoded(GtReads2Twobit *r2t,
    GtUword seqnum, GtTwobitencoding *outputbuffer,
    GtTwobitencoding outputoffset, GtTwobitencoding *lastcodeoffsetptr)
{
  GtUword firstpos, firstcodeidx, lastpos, lastcodeidx, seqlen;
  GtTwobitencoding  inputoffset;

  firstpos = (seqnum == 0)
    ? 0
    : (r2t->seqlen_eqlen > 0)
        ? r2t->seqlen_eqlen * seqnum
        : r2t->seppos[seqnum - 1UL] + 1UL;
  firstcodeidx = GT_DIVBYUNITSIN2BITENC(firstpos);

  lastpos = (r2t->seqlen_eqlen > 0)
    ? r2t->seqlen_eqlen * (seqnum + 1UL) - 1UL
    : r2t->seppos[seqnum];
  lastcodeidx = GT_DIVBYUNITSIN2BITENC(lastpos);

  seqlen = (r2t->seqlen_eqlen > 0)
    ? r2t->seqlen_eqlen
    : lastpos - firstpos + 1UL;

  inputoffset = (GtTwobitencoding)GT_MODBYUNITSIN2BITENC(firstpos);

  if (inputoffset == outputoffset)
  {
    gt_reads2twobit_write_encoded_nocodesshift(r2t, outputbuffer, inputoffset,
        firstcodeidx, lastcodeidx);
    *lastcodeoffsetptr = (GtTwobitencoding)GT_MODBYUNITSIN2BITENC(lastpos + 1);
  }
  else if (inputoffset > outputoffset)
  {
    gt_reads2twobit_write_encoded_leftcodesshift(r2t,
        outputbuffer, inputoffset, outputoffset, firstcodeidx, lastcodeidx);
    *lastcodeoffsetptr = GT_MODBYUNITSIN2BITENC(outputoffset + seqlen);
  }
  else
  {
    gt_reads2twobit_write_encoded_rightcodesshift(r2t,
        outputbuffer, inputoffset, outputoffset, firstcodeidx, lastcodeidx);
    *lastcodeoffsetptr = GT_MODBYUNITSIN2BITENC(outputoffset + seqlen);
  }
  return outputbuffer + GT_DIVBYUNITSIN2BITENC(outputoffset + seqlen);
}

GtUword gt_reads2twobit_mark_mates_of_contained(GtReads2Twobit *r2t,
    GtBitsequence *list)
{
  GtUword libnum, noflibs = gt_array_size(r2t->collection), nofmarked = 0;
  for (libnum = 0; libnum < noflibs; libnum++)
  {
    GtReadsLibraryInfo *rli = gt_array_get(r2t->collection, libnum);
    if (rli->paired && rli->nofseqs > 0)
    {
      GtUword seqnum,
                    last_seqnum = rli->first_seqnum + rli->nofseqs - 1UL;
      gt_assert(rli->nofseqs % 2 == 0);
      for (seqnum = rli->first_seqnum; seqnum < last_seqnum; seqnum += 2UL)
      {
        if (GT_ISIBITSET(list, seqnum) && !GT_ISIBITSET(list, seqnum + 1))
        {
          GT_SETIBIT(list, seqnum + 1UL);
          nofmarked++;
        }
        else if (GT_ISIBITSET(list, seqnum + 1UL) &&
            !GT_ISIBITSET(list, seqnum))
        {
          GT_SETIBIT(list, seqnum);
          nofmarked++;
        }
      }
    }
  }
  return nofmarked;
}

void gt_reads2twobit_delete_sequences(GtReads2Twobit *r2t, GtBitsequence *list)
{
  GtTwobitencoding outputoffset = 0, *outputbuffer;
  GtUword libnum, seqnum, output_seqnum,
                noflibs = gt_array_size(r2t->collection),
                deleted_sequences = 0, deleted_chars = 0, output_startpos = 0;
  for (libnum = 0; libnum < noflibs; libnum++)
  {
    GtReadsLibraryInfo *rli = gt_array_get(r2t->collection, libnum);
    if (rli->nofseqs > 0)
    {
      GtUword deleted_sequences_in_lib = 0, deleted_chars_in_lib = 0,
                    last_seqnum = rli->first_seqnum + rli->nofseqs - 1UL;
      for (seqnum = rli->first_seqnum; seqnum <= last_seqnum; seqnum++)
      {
        if (!GT_ISIBITSET(list, seqnum))
        {
          if (deleted_sequences > 0)
          {
            gt_assert(seqnum >= deleted_sequences);
            output_seqnum = seqnum - deleted_sequences;
            output_startpos = (r2t->seqlen_eqlen > 0)
              ? r2t->seqlen_eqlen * output_seqnum
              : ((output_seqnum == 0) ? 0 :
              r2t->seppos[output_seqnum - 1UL] + 1UL);
            outputbuffer = r2t->twobitencoding +
              GT_DIVBYUNITSIN2BITENC(output_startpos);
            outputoffset = (GtTwobitencoding)
              GT_MODBYUNITSIN2BITENC(output_startpos);
            (void)gt_reads2twobit_write_encoded(r2t, seqnum, outputbuffer,
                outputoffset, &outputoffset);
            if (r2t->seqlen_eqlen == 0)
            {
              GtUword seqlen;
              gt_assert(seqnum > 0);
              seqlen = r2t->seppos[seqnum] - r2t->seppos[seqnum - 1UL];
              gt_assert(seqlen <= r2t->seqlen_max);
              gt_assert(output_startpos + seqlen <= r2t->total_seqlength);
              r2t->seppos[output_seqnum] = output_startpos + seqlen - 1UL;
            }
          }
        }
        else
        {
          deleted_chars_in_lib +=
            gt_reads2twobit_subtract_from_chardistri(r2t, seqnum);
          deleted_sequences++;
          deleted_sequences_in_lib++;
        }
      }
      deleted_chars += deleted_chars_in_lib;
      gt_assert(deleted_sequences_in_lib <= rli->nofseqs);
      if (rli->paired == true)
        gt_assert(deleted_sequences_in_lib % 2 == 0);
      rli->nofseqs -= deleted_sequences_in_lib;
      gt_assert(deleted_chars_in_lib <= rli->total_seqlength);
      rli->total_seqlength -= deleted_chars_in_lib;
    }
  }
  gt_assert(deleted_sequences <= r2t->nofseqs);
  r2t->nofseqs -= deleted_sequences;
  if (r2t->nofseqs == 0)
  {
    r2t->total_seqlength = 0;
  }
  else
  {
    gt_assert(deleted_chars <= r2t->total_seqlength);
    r2t->total_seqlength -= deleted_chars;
  }
  if (deleted_sequences > 0)
  {
    r2t->twobitencoding = gt_realloc(r2t->twobitencoding,
        sizeof (*r2t->twobitencoding) *
        (GT_DIVBYUNITSIN2BITENC(r2t->total_seqlength) + 2UL));
  }
}

static void gt_reads2twobit_set_separators_to_less_frequent_char(
    GtReads2Twobit *r2t)
{
  GtUword seqnum, pos, codenum, posincode;
  GtTwobitencoding sepcode, code, mask, shift;
  sepcode = gt_reads2twobit_less_frequent_char(r2t);
  if (sepcode != r2t->current_sepcode && r2t->nofseqs > 1UL)
  {
    GtUword from, to;
    gt_log_log("changing sepcode from "GT_WU" to "GT_WU"",
                (GtUword) r2t->current_sepcode,
                (GtUword) sepcode);
    from = r2t->seqlen_eqlen > 0 ? 1UL : 0;
    to = r2t->seqlen_eqlen > 0 ? r2t->nofseqs - 1UL : r2t->nofseqs - 2UL;
    for (seqnum = from; seqnum <= to; seqnum++)
    {
      pos = r2t->seqlen_eqlen > 0
        ? seqnum * r2t->seqlen_eqlen - 1UL
        : r2t->seppos[seqnum];
      codenum = GT_DIVBYUNITSIN2BITENC(pos);
      posincode = GT_MODBYUNITSIN2BITENC(pos);
      code = r2t->twobitencoding[codenum];
      shift = (GtTwobitencoding)GT_MULT2(GT_UNITSIN2BITENC - 1UL - posincode);
      mask = ~((GtTwobitencoding)(3UL) << shift);
      gt_assert((code & ~mask) >> shift ==
          (GtTwobitencoding)r2t->current_sepcode);
      code = (code & mask) | ((GtTwobitencoding)sepcode << shift);
      r2t->twobitencoding[codenum] = code;
    }
    r2t->current_sepcode = sepcode;
  }
}

int gt_reads2twobit_write_descriptions(GtReads2Twobit *r2t,
    GtBitsequence *skip, GtError *err)
{
  FILE *desfp = NULL, *sdsfp = NULL;
  int had_err = 0;
  bool with_newline = (r2t->descsfp || (r2t->descs && (!r2t->clipdes)));
  GtUword i, startpos = 0;
  gt_assert(r2t);
  gt_error_check(err);
  gt_assert(r2t->descs != NULL || r2t->descsfp != NULL);
  desfp = gt_fa_fopen_with_suffix(gt_str_get(r2t->indexname),
      GT_DESTABFILESUFFIX, "wb", err);
  if (desfp == NULL)
      had_err = -1;
  if (!had_err) {
    sdsfp = gt_fa_fopen_with_suffix(gt_str_get(r2t->indexname),
        GT_SDSTABFILESUFFIX, "wb", err);
    if (sdsfp == NULL)
        had_err = -1;
  }
  if (!had_err) {
    char *desc = NULL;
    GtUword len = 0, longestdesc = 0, fin = ~0UL, posbuf;
    if (r2t->descsfp != NULL)
    {
      rewind(r2t->descsfp);
      desc = gt_malloc((size_t)r2t->longestdesc+(size_t)1);
    }
    for (i = 0; i < r2t->n_descs; i++)
    {
      if (r2t->descs != NULL)
      {
        desc = (char*) gt_desc_buffer_get_next(r2t->descs);
      }
      else
      {
        (void)gt_xfgets(desc, (int)r2t->longestdesc+(size_t)1, r2t->descsfp);
      }
      if (skip && GT_ISIBITSET(skip, i))
        continue;
      /* note: do not use longestdesc = r2t->longestdesc (for diskbased)
         or the value stored in the desc_buffer (for membased), as these
         values are not up to date when <skip> is applied */
      len = (GtUword)strlen(desc);
      if (with_newline)
      {
        gt_assert(len > 0);
        len--;
      }
      if (len > longestdesc)
        longestdesc = len;
      if (startpos > 0)
      {
        posbuf = startpos - (GtUword)1;
        gt_xfwrite_one(&posbuf, sdsfp);
      }
      gt_xfputs(desc, desfp);
      if (!with_newline)
        gt_xfputc((int) '\n', desfp);
      startpos += (len+(GtUword)1);
    }
    gt_xfwrite_one(&longestdesc, desfp);
    gt_xfwrite_one(&fin, desfp);
    if (r2t->descsfp != NULL)
    {
      gt_free(desc);
    }
  }
  gt_fa_fclose(desfp);
  gt_fa_fclose(sdsfp);
  return had_err;
}

static int gt_reads2twobit_write_encseq_eqlen(GtReads2Twobit *r2t,
    GtError *err)
{
  int had_err = 0;
  GtFilelengthvalues *filelengthtab;
  GtStrArray *filenametab;

  gt_assert(r2t->seqlen_eqlen > 0);
  gt_reads2twobit_collect_fileinfo(r2t, &filelengthtab, &filenametab);
  gt_reads2twobit_set_separators_to_less_frequent_char(r2t);
  gt_reads2twobit_zeropad_tbe(r2t);
  had_err = gt_encseq_equallength_write_twobitencoding_to_file(
      gt_str_get(r2t->indexname), r2t->total_seqlength,
      r2t->seqlen_eqlen - 1UL, r2t->twobitencoding, r2t->nofseqs,
      gt_array_size(r2t->collection), filelengthtab, filenametab,
      r2t->chardistri, err);
  gt_free(filelengthtab);
  gt_str_array_delete(filenametab);
  return had_err;
}

#define GT_READS2TWOBIT_SIZEOFREP(SAT)\
  gt_encseq_determine_size(SAT, r2t->total_seqlength, r2t->nofseqs,\
      numofdbfiles, lengthofdbfilenames, 0, 4U, 2U, 0);

static int gt_reads2twobit_write_encseq_varlen(GtReads2Twobit *r2t,
    GtError *err)
{
  int had_err = 0;
  GtFilelengthvalues *filelengthtab;
  GtStrArray *filenametab;
  uint64_t sizeofrep, minsizeofrep;
  GtUword numofdbfiles = gt_array_size(r2t->collection);
  GtEncseqAccessType sat;
  GtUword idx, lengthofdbfilenames = 0;

  gt_assert(r2t->seppos != NULL);
  gt_reads2twobit_collect_fileinfo(r2t, &filelengthtab, &filenametab);
  for (idx = 0; idx < gt_str_array_size(filenametab); idx++)
  {
    lengthofdbfilenames += gt_str_length(gt_str_array_get_str(
          filenametab, idx)) + 1UL;
  }
  sat = GT_ACCESS_TYPE_UCHARTABLES;
  minsizeofrep = GT_READS2TWOBIT_SIZEOFREP(sat);
  sizeofrep = GT_READS2TWOBIT_SIZEOFREP(GT_ACCESS_TYPE_USHORTTABLES);
  if (sizeofrep < minsizeofrep)
  {
    sat = GT_ACCESS_TYPE_USHORTTABLES;
    minsizeofrep = sizeofrep;
  }
  sizeofrep = GT_READS2TWOBIT_SIZEOFREP(GT_ACCESS_TYPE_UINT32TABLES);
  if (sizeofrep < minsizeofrep)
    sat = GT_ACCESS_TYPE_UINT32TABLES;
  gt_reads2twobit_set_separators_to_less_frequent_char(r2t);
  gt_reads2twobit_zeropad_tbe(r2t);
  had_err = gt_encseq_seppos2ssptab(gt_str_get(r2t->indexname),
      r2t->total_seqlength, r2t->nofseqs, r2t->seppos, err);
  if (!had_err)
    had_err = gt_encseq_generic_write_twobitencoding_to_file(
        gt_str_get(r2t->indexname), r2t->total_seqlength,
        sat, 0, r2t->seqlen_min - 1UL, r2t->seqlen_max - 1UL, 0, 0,
        r2t->seqlen_max - 1UL, r2t->twobitencoding, r2t->nofseqs,
        numofdbfiles, filelengthtab, filenametab, r2t->chardistri, err);
  gt_free(filelengthtab);
  gt_str_array_delete(filenametab);
  return had_err;
}

int gt_reads2twobit_write_encseq(GtReads2Twobit *r2t, GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  gt_assert(r2t != NULL);
  if (r2t->nofseqs == 0)
  {
    gt_log_log("read set is empty, no encseq was written");
    return had_err;
  }
  gt_assert(r2t->twobitencoding != NULL);
  gt_log_log("write encseq %s", gt_str_get(r2t->indexname));
  gt_log_log("seqlen_eqlen = "GT_WU"", r2t->seqlen_eqlen);
  if (r2t->seqlen_eqlen > 0)
  {
    had_err = gt_reads2twobit_write_encseq_eqlen(r2t, err);
  }
  else
  {
    had_err = gt_reads2twobit_write_encseq_varlen(r2t, err);
  }
  return had_err;
}

void gt_reads2twobit_sort(GtReads2Twobit *r2t, GtCompareWithData cmp,
    void *cmp_data)
{
  GtTwobitencoding *tbe, *tbe_next, offset;
  GtUword i, *order;

  gt_assert(r2t != NULL);

  order = gt_malloc(sizeof (*order) * r2t->nofseqs);
  for (i = 0; i < r2t->nofseqs; i++)
    order[i] = i;

  gt_qsort_r(order, (size_t)r2t->nofseqs, sizeof *order, cmp_data, cmp);

  tbe = gt_malloc(sizeof (*tbe) *
      GT_DIVBYUNITSIN2BITENC(r2t->total_seqlength) + 2UL);

  offset = 0;
  tbe_next = tbe;
  for (i = 0; i < r2t->nofseqs; i++)
  {
    tbe_next = gt_reads2twobit_write_encoded(r2t, order[i], tbe_next, offset,
        &offset);
    if (r2t->seqlen_eqlen == 0)
    {
      /* recycle order memory to store new seppos */
      order[i] = (order[i] == 0)
        ? r2t->seppos[0]
        : r2t->seppos[order[i]] - (r2t->seppos[order[i] - 1UL] + 1UL);
      if (i > 0)
        order[i] += (order[i - 1UL] + 1UL);
    }
  }

  gt_free(r2t->twobitencoding);
  r2t->twobitencoding = tbe;

  if (r2t->seqlen_eqlen == 0)
  {
    gt_free(r2t->seppos);
    r2t->seppos = order;
  }
  else
  {
    gt_free(order);
  }
}

int gt_reads2twobit_write_seppos(GtReads2Twobit *r2t, char* path,
    GtBitsequence *skip, GtError *err)
{
  int had_err = 0;
  FILE *file;
  GtUword pos, seqnum;
  if (r2t->seppos == NULL)
    return 0;
  file = gt_fa_fopen(path, "wb", err);
  if (file == NULL)
    had_err = -1;
  if (!had_err)
  {
    if (!GT_ISIBITSET(skip, 0))
    {
      gt_xfwrite(r2t->seppos, sizeof (r2t->seppos), (size_t)1, file);
      pos = r2t->seppos[0] + 1UL;
    }
    else
      pos = 0;
    for (seqnum = 1UL; seqnum < r2t->nofseqs; seqnum++)
    {
      if (!GT_ISIBITSET(skip, seqnum))
      {
        pos += r2t->seppos[seqnum] - r2t->seppos[seqnum - 1UL] - 1UL;
        gt_xfwrite(&pos, sizeof (pos), (size_t)1, file);
        pos++;
      }
    }
    gt_fa_fclose(file);
  }
  return had_err;
}

GtTwobitencoding *gt_reads2twobit_export_twobitencoding(
    const GtReads2Twobit *r2t)
{
  gt_assert(r2t != NULL);
  return r2t->twobitencoding;
}

GtUword *gt_reads2twobit_export_seppos(const GtReads2Twobit *r2t)
{
  gt_assert(r2t != NULL);
  return r2t->seppos;
}

GtUword gt_reads2twobit_nof_invalid_seqs(const GtReads2Twobit *r2t)
{
  gt_assert(r2t != NULL);
  return r2t->invalid_sequences;
}

GtUword gt_reads2twobit_invalid_seqs_totallength(
    const GtReads2Twobit *r2t)
{
  gt_assert(r2t != NULL);
  return r2t->invalid_total_length;
}

GtUword gt_reads2twobit_nofseqs(const GtReads2Twobit *r2t)
{
  gt_assert(r2t != NULL);
  return r2t->nofseqs;
}

GtUword gt_reads2twobit_seqlen_eqlen(const GtReads2Twobit *r2t)
{
  gt_assert(r2t != NULL);
  return r2t->seqlen_eqlen;
}

GtUword gt_reads2twobit_seqlen_max(const GtReads2Twobit *r2t)
{
  gt_assert(r2t != NULL);
  return r2t->seqlen_max;
}

GtUword gt_reads2twobit_seqlen_min(const GtReads2Twobit *r2t)
{
  gt_assert(r2t != NULL);
  return r2t->seqlen_min;
}

GtUword gt_reads2twobit_total_seqlength(const GtReads2Twobit *r2t)
{
  gt_assert(r2t != NULL);
  return r2t->total_seqlength;
}

void gt_reads2twobit_write_libraries_table(const GtReads2Twobit *r2t,
    FILE *rlt_fp)
{
  GtUword noflibs, lnum;
  GtReadsLibrariesTable *rlt;
  noflibs = gt_array_size(r2t->collection);
  gt_assert(noflibs > 0);
  rlt = gt_reads_libraries_table_new(noflibs);
  for (lnum = 0; lnum < noflibs; lnum++)
  {
    GtReadsLibraryInfo *rli;
    rli = gt_array_get(r2t->collection, lnum);
    gt_reads_libraries_table_add(rlt, rli->first_seqnum, rli->insertlength,
        rli->stdev, rli->paired);
  }
  gt_reads_libraries_table_save(rlt, rlt_fp);
  gt_reads_libraries_table_delete(rlt);
}

void gt_reads2twobit_write_hplengths(const GtReads2Twobit *r2t, FILE *out_fp)
{
  gt_assert(r2t != NULL);
  gt_assert(r2t->hplengths != NULL);
  gt_assert(out_fp != NULL);
  gt_hplstore_save(r2t->hplengths, out_fp);
}

/* approx because I use also homopolymers in skipped sequences */
double gt_reads2twobit_approx_average_hplength(const GtReads2Twobit *r2t)
{
  return r2t->approx_avhlen;
}

void gt_reads2twobit_enable_descs(GtReads2Twobit *r2t, bool clipped,
                                  bool membased)
{
  gt_assert(r2t != NULL);
  r2t->clipdes = clipped;
  if (membased)
  {
    r2t->descs = gt_desc_buffer_new();
    if (clipped)
      gt_desc_buffer_set_clip_at_whitespace(r2t->descs);
  }
  else
  {
    r2t->descsfp = gt_xtmpfp_generic(NULL, TMPFP_AUTOREMOVE);
  }
}
