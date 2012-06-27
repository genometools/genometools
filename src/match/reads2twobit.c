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

#ifndef S_SPLINT_S
#include <ctype.h>
#endif
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
#include "core/splitter_api.h"
#include "core/xansi_api.h"
#include "match/reads_library.h"
#include "match/reads2twobit.h"

#define GT_READS2TWOBIT_ALPHASIZE 4U

typedef struct {
  GtStr *filename1;
  GtStr *filename2; /* NULL for unpaired */
  unsigned long insertlength;
  unsigned long total_filelength;
  unsigned long total_seqlength;
  unsigned long first_seqnum;
  unsigned long nofseqs;
} GtReadsLibraryInfo;

struct GtReads2Twobit
{
  GtStr *indexname;
  GtArray *collection;
  GtTwobitencoding *twobitencoding;
  unsigned long chardistri[GT_READS2TWOBIT_ALPHASIZE];
  unsigned long *seppos;
  unsigned long nofseqs;
  unsigned long seqlen_eqlen, seqlen_max, seqlen_min;
  unsigned long total_seqlength;
  GtTwobitencoding current_sepcode;
  unsigned long invalid_sequences, invalid_total_length;
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
  return r2t;
}

void gt_reads2twobit_delete(GtReads2Twobit *r2t)
{
  if (r2t != NULL)
  {
    unsigned long lnum, noflibs;
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
    gt_free(r2t->twobitencoding);
    gt_free(r2t->seppos);
    gt_free(r2t);
  }
}

void gt_reads2twobit_add_paired(GtReads2Twobit *r2t, const GtStr *filename1,
    const GtStr *filename2, unsigned long insertlength)
{
  GtReadsLibraryInfo rli;
  gt_assert(r2t != NULL);
  gt_assert(filename1 != NULL);
  rli.filename1 = gt_str_clone(filename1);
  if (filename2 != NULL)
    rli.filename2 = gt_str_clone(filename2);
  else
    rli.filename2 = NULL;
  rli.insertlength = insertlength;
  rli.total_seqlength = 0;
  rli.first_seqnum = 0;
  rli.nofseqs = 0;
  rli.total_filelength = (unsigned long)gt_file_size(gt_str_get(filename1));
  if (filename2 != NULL)
    rli.total_filelength += (unsigned long)gt_file_size(gt_str_get(filename2));
  gt_array_add(r2t->collection, rli);
}

int gt_reads2twobit_add_library(GtReads2Twobit *r2t, const GtStr *libspec,
                                GtError *err)
{
  GtSplitter *s;
  char *libspec_copy;
  int had_err = 0;

  gt_assert(r2t != NULL);
  gt_assert(libspec != NULL);
  gt_assert(gt_str_length(libspec) > 0);
  libspec_copy = gt_malloc(sizeof (*libspec_copy) *
      (gt_str_length(libspec) + 1UL));
  (void)strcpy(libspec_copy, gt_str_get(libspec));
  s = gt_splitter_new();
  gt_splitter_split(s, libspec_copy, gt_str_length(libspec),
      GT_READS2TWOBIT_LIBNAMESEP);
  gt_log_log("reads2twobit add library: %s", gt_str_get(libspec));
  if (gt_splitter_size(s) == 1UL)
  {
    gt_reads2twobit_add_unpaired(r2t, libspec);
  }
  else if (gt_splitter_size(s) == 3UL)
  {
    GtStr *filename1 = gt_str_new_cstr(gt_splitter_get_token(s, 0)),
          *filename2 = gt_str_new_cstr(gt_splitter_get_token(s, 1UL));
    const char *insertlength_str = gt_splitter_get_token(s, 2UL);
    unsigned long insertlength;
    had_err = gt_parse_ulong(&insertlength, insertlength_str);
    if (had_err != 0)
    {
      gt_error_set(err, "insertlength specification not valid: %s",
          insertlength_str);
    }
    if (!had_err)
    {
      gt_reads2twobit_add_paired(r2t, filename1, filename2, insertlength);
    }
    gt_str_delete(filename1);
    gt_str_delete(filename2);
  }
  else
  {
    gt_error_set(err, "library specification not valid: %s\n"
        "use a single filename or filename1%cfilename2%cinsertlength\n"
        "filenames containing %c are not allowed", gt_str_get(libspec),
      GT_READS2TWOBIT_LIBNAMESEP, GT_READS2TWOBIT_LIBNAMESEP,
      GT_READS2TWOBIT_LIBNAMESEP);
    had_err = -1;
  }
  gt_free(libspec_copy);
  gt_splitter_delete(s);
  return had_err;
}

#define GT_READS2TWOBIT_CODE_UNDEF ((GtTwobitencoding)ULONG_MAX)

typedef struct {
  GtTwobitencoding *tbe_next;
  GtTwobitencoding kmercode;
  unsigned short codepos;
  unsigned long chardistri[GT_READS2TWOBIT_ALPHASIZE];
  unsigned long globalpos;
  unsigned long nofseqs;
  unsigned long seppos_nextfree;
  unsigned long seqlen_max;
  unsigned long seqlen_min;
} GtReads2TwobitEncodeInfo;

#define GT_READS2TWOBIT_INIT_ENCODE_INFO(EI, TBE)\
  memset(&(EI), 0, sizeof (EI));\
  (EI).tbe_next = (TBE)

#define GT_READS2TWOBIT_COPY_ENCODE_INFO(EISRC, EIDEST)\
  memcpy(&(EIDEST), &(EISRC), sizeof (EISRC))

typedef struct {
  GtTwobitencoding char2code[UCHAR_MAX + 1];
  unsigned long inputfiles_totallength;
  unsigned long seqlen, seqlen_first;
  GtReads2TwobitEncodeInfo current, backup;
  bool varlen_mode, invalid_mode;
  unsigned long invalid_sequences, invalid_total_length;
  unsigned long *seppos;
  unsigned long seppos_alloc;
} GtReads2TwobitEncodeState;

static void gt_reads2twobit_init_encode(GtReads2Twobit *r2t,
    GtReads2TwobitEncodeState *state)
{
  const unsigned long noflibraries = gt_array_size(r2t->collection);
  unsigned long i;
  state->inputfiles_totallength = 0;
  for (i = 0; i < noflibraries; i++)
  {
    const GtReadsLibraryInfo *rli = gt_array_get(r2t->collection, i);
    state->inputfiles_totallength += rli->total_filelength;
  }
  r2t->twobitencoding = gt_malloc(sizeof (*r2t->twobitencoding) *
    GT_DIVBYUNITSIN2BITENC(state->inputfiles_totallength) + 2UL);
  GT_READS2TWOBIT_INIT_ENCODE_INFO(state->current, r2t->twobitencoding);
  GT_READS2TWOBIT_COPY_ENCODE_INFO(state->current, state->backup);
  for (i = 0; i < (unsigned long)(UCHAR_MAX) + 1UL; i++)
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
  state->seqlen_first = 0;
  state->varlen_mode = false;
  state->invalid_mode = false;
  state->invalid_sequences = 0;
  state->invalid_total_length = 0;
  gt_assert(r2t->seppos == NULL);
  state->seppos = NULL;
}

#define GT_READS2TWOBIT_READBUFFER_SIZE ((size_t)256)

#define GT_READS2TWOBIT_WRITECODE_NOCOUNT(EI, CODE, LEN) \
  (LEN)++;\
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
    unsigned long currentpos)
{
  gt_assert(state->seppos == NULL);
  gt_assert(currentpos > 0);
  gt_assert(state->inputfiles_totallength > currentpos);
  state->seppos_alloc = state->inputfiles_totallength / currentpos;
  gt_log_log("rough estimate of nofseqs = %lu", state->seppos_alloc);
  state->seppos = gt_malloc(sizeof (state->seppos) * state->seppos_alloc);
  state->current.seppos_nextfree = 0;
}

#define GT_READS2TWOBIT_SEPPOS_INC ((size_t)(1 << 14))

static void gt_reads2twobit_append_seppos(GtReads2TwobitEncodeState *state)
{
  gt_assert(state->seppos != NULL);
  if (state->current.seppos_nextfree == state->seppos_alloc)
  {
    (state->seppos_alloc) += GT_READS2TWOBIT_SEPPOS_INC;
    state->seppos = gt_realloc(state->seppos,
        sizeof (state->seppos) * state->seppos_alloc);
  }
  gt_assert(state->current.globalpos > 0);
  state->seppos[state->current.seppos_nextfree] =
    state->current.globalpos - 1UL;
  state->current.seppos_nextfree++;
}

static void gt_reads2twobit_switch_to_varlen_mode(
    GtReads2TwobitEncodeState *state)
{
  unsigned long seqnum;
  gt_assert(state->varlen_mode == false);
  state->varlen_mode = true;
  gt_assert(state->current.nofseqs > 1UL);
  gt_assert(state->seqlen_first != state->seqlen);
  gt_assert(state->seqlen > 1UL);
  gt_log_log("readset is varlen: sequences 0..%lu are "
      "%lu bp long, sequence %lu is %lu bp long",
      state->current.nofseqs - 2UL, state->seqlen_first - 1UL,
      state->current.nofseqs - 1UL, state->seqlen - 1UL);
  gt_assert(state->current.globalpos == 0);
  gt_reads2twobit_init_seppos(state,
      state->seqlen_first * (state->current.nofseqs - 2UL) +
      state->seqlen);
  for (seqnum = 0; seqnum < state->current.nofseqs - 1UL; seqnum++)
  {
    state->current.globalpos += state->seqlen_first;
    gt_reads2twobit_append_seppos(state);
  }
  state->current.globalpos += state->seqlen;
  gt_reads2twobit_append_seppos(state);
  gt_assert(state->current.seppos_nextfree == state->current.nofseqs);
  state->current.seqlen_max = MAX(state->seqlen_first, state->seqlen);
  state->current.seqlen_min = MIN(state->seqlen_first, state->seqlen);
  state->seqlen_first = 0;
}

#define GT_READS2TWOBIT_DEFAULT_SEPARATOR (GtTwobitencoding)3

static inline void gt_reads2twobit_process_sequence_end(
    GtReads2TwobitEncodeState *state)
{
  GT_READS2TWOBIT_WRITECODE_NOCOUNT(state->current,
      GT_READS2TWOBIT_DEFAULT_SEPARATOR, state->seqlen);
  gt_assert(!state->invalid_mode);
  if (state->varlen_mode)
  {
    if (state->seqlen > state->current.seqlen_max)
      state->current.seqlen_max = state->seqlen;
    if (state->seqlen < state->current.seqlen_min)
      state->current.seqlen_min = state->seqlen;
    state->current.globalpos += state->seqlen;
    gt_reads2twobit_append_seppos(state);
  }
  else
  {
    if (state->current.nofseqs > 1UL)
    {
      if (state->seqlen != state->seqlen_first)
        gt_reads2twobit_switch_to_varlen_mode(state);
    }
    else
      state->seqlen_first = state->seqlen;
  }
}

static inline void gt_reads2twobit_prepare_for_new_sequence(
    GtReads2TwobitEncodeState *state)
{
  GT_READS2TWOBIT_COPY_ENCODE_INFO(state->current, state->backup);
  state->current.nofseqs++;
  state->seqlen = 0;
  state->invalid_mode = false;
}

#define GT_READS2TWOBIT_SKIP_TO_DESCRIPTION_END(LINE, FILEPTR, FGETSRETVAL)\
  /* handle the case in which a description is longer than the line buffer: */\
  while (strlen(LINE) == GT_READS2TWOBIT_READBUFFER_SIZE - (size_t) 1 \
      && (LINE)[GT_READS2TWOBIT_READBUFFER_SIZE - 2] != '\n' \
      && (FGETSRETVAL) != NULL) \
  {\
    (FGETSRETVAL) = \
      fgets((LINE), (int)GT_READS2TWOBIT_READBUFFER_SIZE, (FILEPTR));\
  }

static void gt_reads2twobit_switch_to_invalid_mode(
    GtReads2TwobitEncodeState *state)
{
  state->invalid_mode = true;
  state->invalid_sequences++;
  state->invalid_total_length += state->seqlen;
  GT_READS2TWOBIT_COPY_ENCODE_INFO(state->backup, state->current);
}

static inline void gt_reads2twobit_process_sequence_line(
    GtReads2TwobitEncodeState *state, const char *line)
{
  unsigned long j = 0;
  GtTwobitencoding nextcode;
  char c;
  while (true)
  {
    c = line[j++];
    if (!state->invalid_mode && (nextcode = state->char2code[(unsigned char)c])
        != GT_READS2TWOBIT_CODE_UNDEF)
    {
      GT_READS2TWOBIT_WRITECODE(state->current, nextcode, state->seqlen);
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
    const unsigned long invalid_tl_before = state->invalid_total_length,
          invalid_s_before = state->invalid_sequences;
    char line[GT_READS2TWOBIT_READBUFFER_SIZE], *fgetsretval;
    rli->first_seqnum = state->current.nofseqs;
    while ((fgetsretval = fgets(line, (int)GT_READS2TWOBIT_READBUFFER_SIZE,
              file)) == line)
    {
      if (line[0] == '>')
      {
        if (state->current.nofseqs > rli->first_seqnum && !state->invalid_mode)
        {
          gt_reads2twobit_process_sequence_end(state);
        }
        GT_READS2TWOBIT_SKIP_TO_DESCRIPTION_END(line, file, fgetsretval);
        gt_reads2twobit_prepare_for_new_sequence(state);
      }
      else if (!state->invalid_mode)
      {
        gt_reads2twobit_process_sequence_line(state, line);
      }
    }
    if (!state->invalid_mode)
      gt_reads2twobit_process_sequence_end(state);
    gt_assert(state->current.nofseqs >= rli->first_seqnum);
    rli->nofseqs = state->current.nofseqs - rli->first_seqnum;
    rli->total_seqlength = state->varlen_mode
      ?  state->seppos[state->current.nofseqs - 1UL] + 1UL -
         ((rli->first_seqnum == 0) ? 0 : state->seppos[rli->first_seqnum - 1UL])
      : state->seqlen_first * rli->nofseqs;
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

static inline void gt_reads2twobit_process_mate_pair(
    GtReads2TwobitEncodeState *state, char *line2, FILE *file2)
{
  char *fgetsretval = NULL;
  bool was_invalid = state->invalid_mode;
  unsigned long prev_seqlen = state->seqlen;
  if (line2[0] != '>')
  {
    fgetsretval = fgets(line2, (int)GT_READS2TWOBIT_READBUFFER_SIZE,
        file2);
    gt_assert(fgetsretval == line2);
  }
  gt_assert(line2[0] == '>');
  GT_READS2TWOBIT_SKIP_TO_DESCRIPTION_END(line2, file2, fgetsretval);
  if (!state->invalid_mode)
  {
    state->current.nofseqs++;
    state->seqlen = 0;
  }
  else
  {
    state->invalid_sequences++;
  }
  while (true)
  {
    fgetsretval = fgets(line2,
          (int)GT_READS2TWOBIT_READBUFFER_SIZE, file2);
    if (fgetsretval != line2 || line2[0] == '>')
      break;
    gt_reads2twobit_process_sequence_line(state, line2);
  }
  if (!state->invalid_mode)
    gt_reads2twobit_process_sequence_end(state);
  if (!was_invalid && state->invalid_mode)
  {
    state->invalid_sequences++;
    state->invalid_total_length += (prev_seqlen - 1UL);
  }
}

static int gt_reads2twobit_encode_paired_library(
    GtReads2TwobitEncodeState *state, GtReadsLibraryInfo *rli, GtError *err)
{
  int had_err = 0;
  FILE *file1, *file2 = NULL;
  file1 = gt_fa_fopen(gt_str_get(rli->filename1), "r", err);
  if (file1 == NULL)
    had_err = -1;
  if (!had_err)
  {
    file2 = gt_fa_fopen(gt_str_get(rli->filename2), "r", err);
    if (file2 == NULL)
      had_err = -1;
  }
  if (!had_err)
  {
    const unsigned long invalid_tl_before = state->invalid_total_length,
          invalid_s_before = state->invalid_sequences;
    char line1[GT_READS2TWOBIT_READBUFFER_SIZE],
         line2[GT_READS2TWOBIT_READBUFFER_SIZE], *fgetsretval;
    rli->first_seqnum = state->current.nofseqs;
    while ((fgetsretval = fgets(line1, (int)GT_READS2TWOBIT_READBUFFER_SIZE,
              file1)) == line1)
    {
      if (line1[0] == '>')
      {
        if (state->current.nofseqs > rli->first_seqnum)
        {
          if (!state->invalid_mode)
            gt_reads2twobit_process_sequence_end(state);
          gt_reads2twobit_process_mate_pair(state, line2, file2);
        }
        else
        {
          if (state->invalid_mode)
            gt_reads2twobit_process_mate_pair(state, line2, file2);
        }
        GT_READS2TWOBIT_SKIP_TO_DESCRIPTION_END(line1, file1, fgetsretval);
        gt_reads2twobit_prepare_for_new_sequence(state);
      }
      else if (!state->invalid_mode)
      {
        gt_reads2twobit_process_sequence_line(state, line1);
      }
    }
    if (!state->invalid_mode)
      gt_reads2twobit_process_sequence_end(state);
    gt_reads2twobit_process_mate_pair(state, line2, file2);
    gt_assert(state->current.nofseqs >= rli->first_seqnum);
    rli->nofseqs = state->current.nofseqs - rli->first_seqnum;
    rli->total_seqlength = state->varlen_mode
      ?  state->seppos[state->current.nofseqs - 1UL] + 1UL -
         ((rli->first_seqnum == 0) ? 0 : state->seppos[rli->first_seqnum - 1UL])
      : state->seqlen_first * rli->nofseqs;
    /* the following is not necessary, but is useful for the tests */
    rli->total_filelength -= (state->invalid_total_length - invalid_tl_before +
        3UL * (state->invalid_sequences - invalid_s_before));
  }
  if (!had_err)
    had_err = gt_reads2twobit_close_file(file1, rli->filename1, err);
  else if (file1 != NULL)
    gt_fa_fclose(file1);
  if (!had_err)
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
    unsigned long shift = (GT_UNITSIN2BITENC - state->current.codepos) << 1UL;
    *(state->current.tbe_next++) =
      (GtTwobitencoding)(state->current.kmercode << shift);
  }
  if (state->current.nofseqs > 0)
  {
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
    r2t->seqlen_eqlen = state->seqlen_first;
    r2t->seqlen_max = state->seqlen_first;
    r2t->seqlen_min = state->seqlen_first;
    r2t->total_seqlength = state->seqlen_first * state->current.nofseqs - 1UL;
    gt_assert(state->seppos == NULL);
  }
  gt_reads2twobit_tbe_flush_and_realloc(r2t, state);
}

int gt_reads2twobit_encode(GtReads2Twobit *r2t, GtError *err)
{
  int had_err = 0;
  const unsigned long noflibraries = gt_array_size(r2t->collection);
  unsigned long libnum;
  GtReads2TwobitEncodeState state;

  gt_error_check(err);
  gt_assert(r2t != NULL);
  gt_assert(r2t->twobitencoding == NULL);
  gt_reads2twobit_init_encode(r2t, &state);
  for (libnum = 0; libnum < noflibraries && !had_err; libnum++)
  {
    GtReadsLibraryInfo *rli = gt_array_get(r2t->collection, libnum);
    if (rli->filename2 == NULL)
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
  unsigned long i, noflibraries;
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
      (*filelengthtab)[i].length = (uint64_t)rli->total_filelength;
    }
    libname = gt_str_clone(rli->filename1);
    if (rli->filename2 != NULL)
    {
      gt_str_append_char(libname, GT_READS2TWOBIT_LIBNAMESEP);
      gt_str_append_str(libname, rli->filename2);
      gt_str_append_char(libname, GT_READS2TWOBIT_LIBNAMESEP);
      gt_str_append_ulong(libname, rli->insertlength);
    }
    gt_str_array_add(*filenametab, libname);
    gt_str_delete(libname);
  }
}

static inline GtTwobitencoding gt_reads2twobit_less_frequent_char(
    GtReads2Twobit *r2t)
{
  GtTwobitencoding i, code;
  unsigned long lowest_value, value;
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
  gt_log_log("less frequent char code: %lu", (unsigned long)code);
  return code;
}

static void gt_reads2twobit_zeropad_tbe(GtReads2Twobit *r2t)
{
  unsigned long pos, codenum, posincode, shift;
  gt_assert(r2t->total_seqlength > 0);
  pos = r2t->total_seqlength - 1UL;
  codenum = GT_DIVBYUNITSIN2BITENC(pos);
  posincode = GT_MODBYUNITSIN2BITENC(pos);
  if (posincode < (unsigned long)GT_UNITSIN2BITENC - 1UL)
  {
    shift = GT_MULT2(GT_UNITSIN2BITENC - 1UL - posincode);
    r2t->twobitencoding[codenum] =
      (GtTwobitencoding)((r2t->twobitencoding[codenum] >> shift) << shift);
  }
  r2t->twobitencoding[codenum + 1UL] = 0;
}

static void gt_reads2twobit_seek_sequence(const GtReads2Twobit *r2t,
    unsigned long seqnum, unsigned long *seqlen, GtTwobitencoding *firstcode,
    unsigned long *charsinfirstcode, GtTwobitencoding **nextcode_ptr)
{
  unsigned long pos;
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

void gt_reads2twobit_decode_sequence(const GtReads2Twobit *r2t,
    unsigned long seqnum, char *decoded)
{
  GtTwobitencoding code;
  unsigned long pos, seqlen, charsincode;
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
      charsincode = (unsigned long)GT_UNITSIN2BITENC;
    }
    *(nextdecoded++) = code2char[code >> ((--charsincode) << 1) & 3];
  }
  *(nextdecoded++) = '\n';
  *(nextdecoded) = '\0';
}

static unsigned long gt_reads2twobit_subtract_from_chardistri(
    GtReads2Twobit *r2t, unsigned long seqnum)
{
  GtTwobitencoding code;
  unsigned long pos, seqlen, charsincode;
  GtTwobitencoding *nextencoded;

  gt_reads2twobit_seek_sequence(r2t, seqnum, &seqlen, &code, &charsincode,
      &nextencoded);
  for (pos = 0; pos < seqlen - 1UL; pos++)
  {
    if (charsincode == 0)
    {
      code = *(nextencoded++);
      charsincode = (unsigned long)GT_UNITSIN2BITENC;
    }
    r2t->chardistri[code >> ((--charsincode) << 1) & 3]--;
  }
  return seqlen;
}

void gt_reads2twobit_decode_range(const GtReads2Twobit *r2t,
    GtFile *outfp, unsigned long seqnum_from, unsigned long nofseqs,
    const GtBitsequence *skip)
{
  GtTwobitencoding code;
  unsigned short charsincode;
  const char code2char[] = "acgt";
  unsigned long seqnum, pos, nextsep, nextdecoded, seqnum_to;
  const GtTwobitencoding *nextencoded = r2t->twobitencoding;
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
    unsigned long firstcodeidx, unsigned long lastcodeidx)
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
    unsigned long firstcodeidx, unsigned long lastcodeidx)
{
  const GtTwobitencoding netoffset = inputoffset - outputoffset;
  const GtTwobitencoding shiftright =
    GT_MULT2(GT_UNITSIN2BITENC - netoffset);
  const GtTwobitencoding shiftleft = GT_MULT2(netoffset);
  unsigned long i;
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
    unsigned long firstcodeidx, unsigned long lastcodeidx)
{
  const GtTwobitencoding netoffset = outputoffset - inputoffset;
  const GtTwobitencoding shiftright = GT_MULT2(netoffset);
  const GtTwobitencoding shiftleft =
    GT_MULT2(GT_UNITSIN2BITENC - netoffset);
  unsigned long i;
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

unsigned long gt_reads2twobit_write_encoded(GtReads2Twobit *r2t,
    unsigned long seqnum, GtTwobitencoding *outputbuffer,
    GtTwobitencoding outputoffset, GtTwobitencoding *lastcodeoffsetptr)
{
  unsigned long firstpos, firstcodeidx, lastpos, lastcodeidx, seqlen;
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
  return GT_DIVBYUNITSIN2BITENC(seqlen + outputoffset + GT_UNITSIN2BITENC - 1);
}

static inline void gt_reads2twobit_handle_deleted_mates(
    GT_UNUSED GtReads2Twobit *r2t, GtReadsLibraryInfo *rli, GtBitsequence *list)
{
  unsigned long seqnum, last_seqnum = rli->first_seqnum + rli->nofseqs - 1UL;
  gt_assert(rli->nofseqs % 2 == 0);
  gt_assert(rli->nofseqs > 0);
  for (seqnum = rli->first_seqnum; seqnum < last_seqnum; seqnum += 2UL)
  {
    if (GT_ISIBITSET(list, seqnum))
    {
      GT_SETIBIT(list, seqnum + 1UL);
    }
    else if (GT_ISIBITSET(list, seqnum + 1UL))
    {
      GT_SETIBIT(list, seqnum);
    }
  }
}

void gt_reads2twobit_delete_sequences(GtReads2Twobit *r2t, GtBitsequence *list)
{
  GtTwobitencoding outputoffset = 0, *outputbuffer;
  unsigned long libnum, seqnum, output_seqnum,
                noflibs = gt_array_size(r2t->collection),
                deleted_sequences = 0, deleted_chars = 0, output_startpos = 0;
  for (libnum = 0; libnum < noflibs; libnum++)
  {
    GtReadsLibraryInfo *rli = gt_array_get(r2t->collection, libnum);
    if (rli->nofseqs > 0)
    {
      unsigned long deleted_sequences_in_lib = 0, deleted_chars_in_lib = 0,
                    last_seqnum = rli->first_seqnum + rli->nofseqs - 1UL;
      if (rli->filename2 != NULL)
        gt_reads2twobit_handle_deleted_mates(r2t, rli, list);
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
              unsigned long seqlen;
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
      if (rli->filename2 != NULL)
        gt_assert(deleted_sequences_in_lib % 2 == 0);
      rli->nofseqs -= deleted_sequences_in_lib;
      gt_assert(deleted_chars_in_lib <= rli->total_seqlength);
      rli->total_seqlength -= deleted_chars_in_lib;
    }
  }
  gt_assert(deleted_sequences <= r2t->nofseqs);
  r2t->nofseqs -= deleted_sequences;
  gt_assert(deleted_chars <= r2t->total_seqlength);
  r2t->total_seqlength -= deleted_chars;
  if (deleted_sequences > 0)
  {
    r2t->twobitencoding = gt_realloc(r2t->twobitencoding,
        sizeof (*r2t->twobitencoding) *
        (GT_DIVBYUNITSIN2BITENC(r2t->total_seqlength) + 2UL));
  }
}

static void gt_reads2twobit_eqlen_set_separators_to_less_frequent_char(
    GtReads2Twobit *r2t)
{
  unsigned long seqnum, pos, codenum, posincode;
  GtTwobitencoding sepcode, code, mask, shift;
  gt_assert(r2t->seqlen_eqlen > 0);
  sepcode = gt_reads2twobit_less_frequent_char(r2t);
  if (sepcode != r2t->current_sepcode)
  {
    gt_log_log("changing sepcode from %lu to %lu", r2t->current_sepcode,
        sepcode);
    for (seqnum = 1UL; seqnum < r2t->nofseqs; seqnum++)
    {
      pos = seqnum * r2t->seqlen_eqlen - 1UL;
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

static int gt_reads2twobit_write_encseq_eqlen(GtReads2Twobit *r2t,
    GtError *err)
{
  int had_err = 0;
  GtFilelengthvalues *filelengthtab;
  GtStrArray *filenametab;

  gt_assert(r2t->seqlen_eqlen > 0);
  gt_reads2twobit_collect_fileinfo(r2t, &filelengthtab, &filenametab);
  gt_reads2twobit_eqlen_set_separators_to_less_frequent_char(r2t);
  gt_reads2twobit_zeropad_tbe(r2t);
  had_err = gt_encseq_equallength_write_twobitencoding_to_file(
      gt_str_get(r2t->indexname),
      r2t->total_seqlength, r2t->seqlen_eqlen - 1UL, r2t->twobitencoding,
      r2t->nofseqs, gt_array_size(r2t->collection), filelengthtab,
      filenametab, r2t->chardistri, err);
  gt_free(filelengthtab);
  gt_str_array_delete(filenametab);
  return had_err;
}

static int gt_reads2twobit_write_encseq_varlen(const GtReads2Twobit *r2t,
    GtError *err)
{
  GtEncseqEncoder *encoder;
  GtStrArray *infiles;
  GtStrArray *filenametab;
  GtStr *fas_path;
  int had_err = 0;
  unsigned long i, noflibraries = gt_array_size(r2t->collection);

  gt_reads2twobit_collect_fileinfo(r2t, NULL, &filenametab);
  infiles = gt_str_array_new();
  fas_path = gt_str_new();
  for (i = 0; i < noflibraries; i++)
  {
    const GtReadsLibraryInfo *rli = gt_array_get(r2t->collection, i);
    GtFile *fp;
    char *basename;
    gt_str_set(fas_path, "tmpfile.");
    basename = gt_basename(gt_str_array_get(filenametab, i));
    gt_str_append_cstr(fas_path, basename);
    fp = gt_file_new(gt_str_get(fas_path), "w", err);
    if (fp == NULL)
    {
      had_err = -1;
      break;
    }
    gt_log_log("decode library %lu to %s", i, gt_str_get(fas_path));
    gt_reads2twobit_decode_range(r2t, fp, rli->first_seqnum, rli->nofseqs,
        NULL);
    gt_free(basename);
    gt_file_delete(fp);
    gt_str_array_add(infiles, fas_path);
  }
  gt_str_delete(fas_path);
  encoder = gt_encseq_encoder_new();
  gt_encseq_encoder_disable_description_support(encoder);
  gt_encseq_encoder_disable_md5_support(encoder);
  gt_log_log("encode encseq %s", gt_str_get(r2t->indexname));
  had_err = gt_encseq_encoder_encode(encoder, infiles,
      gt_str_get(r2t->indexname), err);
  gt_encseq_encoder_delete(encoder);
  for (i = 0; i < noflibraries; i++)
  {
    gt_xremove(gt_str_array_get(infiles, i));
  }
  gt_str_array_delete(infiles);
  gt_str_array_delete(filenametab);
  return had_err;
}

int gt_reads2twobit_write_encseq(GtReads2Twobit *r2t, GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  gt_assert(r2t != NULL);
  gt_assert(r2t->twobitencoding != NULL);
  gt_log_log("write encseq %s", gt_str_get(r2t->indexname));
  gt_log_log("seqlen_eqlen = %lu", r2t->seqlen_eqlen);
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

int gt_reads2twobit_write_seppos(GtReads2Twobit *r2t, char* path,
    GtBitsequence *skip, GtError *err)
{
  int had_err = 0;
  FILE *file;
  unsigned long pos, seqnum;
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
  gt_assert(r2t->twobitencoding != NULL);
  return r2t->twobitencoding;
}

unsigned long *gt_reads2twobit_export_seppos(const GtReads2Twobit *r2t)
{
  gt_assert(r2t != NULL);
  gt_assert(r2t->twobitencoding != NULL);
  return r2t->seppos;
}

unsigned long gt_reads2twobit_nof_invalid_seqs(const GtReads2Twobit *r2t)
{
  gt_assert(r2t != NULL);
  gt_assert(r2t->twobitencoding != NULL);
  return r2t->invalid_sequences;
}

unsigned long gt_reads2twobit_invalid_seqs_totallength(
    const GtReads2Twobit *r2t)
{
  gt_assert(r2t != NULL);
  gt_assert(r2t->twobitencoding != NULL);
  return r2t->invalid_total_length;
}

unsigned long gt_reads2twobit_nofseqs(const GtReads2Twobit *r2t)
{
  gt_assert(r2t != NULL);
  gt_assert(r2t->twobitencoding != NULL);
  return r2t->nofseqs;
}

unsigned long gt_reads2twobit_seqlen_eqlen(const GtReads2Twobit *r2t)
{
  gt_assert(r2t != NULL);
  gt_assert(r2t->twobitencoding != NULL);
  return r2t->seqlen_eqlen;
}

unsigned long gt_reads2twobit_seqlen_max(const GtReads2Twobit *r2t)
{
  gt_assert(r2t != NULL);
  gt_assert(r2t->twobitencoding != NULL);
  return r2t->seqlen_max;
}

unsigned long gt_reads2twobit_seqlen_min(const GtReads2Twobit *r2t)
{
  gt_assert(r2t != NULL);
  gt_assert(r2t->twobitencoding != NULL);
  return r2t->seqlen_min;
}

unsigned long gt_reads2twobit_total_seqlength(const GtReads2Twobit *r2t)
{
  gt_assert(r2t != NULL);
  gt_assert(r2t->twobitencoding != NULL);
  return r2t->total_seqlength;
}

int gt_reads2twobit_write_libraries_table(const GtReads2Twobit *r2t,
    char *path, GtError *err)
{
  int had_err = 0;
  unsigned long noflibs, lnum;
  GtReadsLibrary *lib_table;
  noflibs = gt_array_size(r2t->collection);
  gt_assert(noflibs > 0);
  lib_table = gt_malloc(sizeof (*lib_table) * noflibs);
  for (lnum = 0; lnum < noflibs; lnum++)
  {
    GtReadsLibraryInfo *rli;
    rli = gt_array_get(r2t->collection, lnum);
    lib_table[lnum].first_seqnum = rli->first_seqnum;
    lib_table[lnum].nofseqs = rli->nofseqs;
    lib_table[lnum].paired = (rli->filename2 != NULL) ? true : false;
    lib_table[lnum].insertlength = rli->insertlength;
  }
  had_err = gt_reads_library_table_write(lib_table, noflibs, path, err);
  gt_free(lib_table);
  return had_err;
}
