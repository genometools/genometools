/*
  Copyright (c) 2010-2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2010-2011 Center for Bioinformatics, University of Hamburg

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

#include <stdio.h>
#include <limits.h>
#include "core/fa.h"
#include "core/fileutils.h"
#include "core/log.h"
#include "core/xansi_api.h"
#include "match/rdj-cntlist.h"

#define GT_CNTLIST_BIT_HEADER   (int)'\0'
#define GT_CNTLIST_BIN_HEADER   (int)'\1'
#define GT_CNTLIST_ASCII_HEADER (int)'['

static inline void gt_cntlist_show_ascii(GtBitsequence *cntlist,
    unsigned long nofreads, FILE *file)
{
  unsigned long i;
  gt_assert(file != NULL);
  fprintf(file, "[n: %lu]\n", nofreads);
  for (i = 0; i < nofreads; i++)
    if (GT_ISIBITSET(cntlist, i))
      fprintf(file, "%lu\n", i);
}

void gt_cntlist_write_bin_header(unsigned long nofreads, FILE *file)
{
  gt_assert(file != NULL);
  gt_xfputc(GT_CNTLIST_BIN_HEADER, file);
  gt_xfputc((char)sizeof(unsigned long), file);
  gt_xfwrite(&(nofreads), sizeof (unsigned long), (size_t)1, file);
}

static inline void gt_cntlist_show_bit(GtBitsequence *cntlist,
    unsigned long nofreads, FILE *file)
{
  gt_assert(file != NULL);
  gt_xfputc(GT_CNTLIST_BIT_HEADER, file);
  gt_xfputc((char)sizeof(unsigned long), file);
  gt_xfwrite(&(nofreads), sizeof (unsigned long), (size_t)1, file);
  gt_xfwrite(cntlist, sizeof (GtBitsequence), GT_NUMOFINTSFORBITS(nofreads),
      file);
}

int gt_cntlist_show(GtBitsequence *cntlist, unsigned long nofreads,
    const char *path, bool binary, GtError *err)
{
  FILE *file;
  gt_assert(cntlist != NULL);
  if (path == NULL)
    file = stdout;
  else
  {
    file = gt_fa_fopen(path, binary ? "wb" : "w", err);
    if (file == NULL)
      return -1;
  }
  gt_assert(file != NULL);
  (binary ? gt_cntlist_show_bit : gt_cntlist_show_ascii)
    (cntlist, nofreads, file);
  if (path != NULL)
    gt_fa_fclose(file);
  return 0;
}

static int gt_cntlist_parse_bin_or_bit_header(FILE *infp,
    unsigned long *nofreads, GtError *err)
{
  int c;
  size_t n;

  gt_assert(infp != NULL && nofreads != NULL);
  gt_error_check(err);
  c = gt_xfgetc(infp);
  if (c == EOF)
  {
    gt_error_set(err, "contained reads list: unexpected end of file");
    return -1;
  }
  else if (c != (char)sizeof(unsigned long))
  {
    gt_error_set(err, "contained reads list: %dbit version "
        "of GenomeTools required to use this list", c * CHAR_BIT);
    return -1;
  }
  n = fread(nofreads, sizeof (unsigned long), (size_t)1, infp);
  if (n != (size_t)1 || *nofreads == 0)
  {
    gt_error_set(err, "contained reads list: unrecognized format");
    return -1;
  }
  return 0;
}

static int gt_cntlist_parse_bit(FILE *infp, bool alloc_cntlist,
    GtBitsequence **cntlist, unsigned long *nofreads, GtError *err)
{
  int had_err = gt_cntlist_parse_bin_or_bit_header(infp, nofreads, err);
  if (had_err == 0)
  {
    size_t n;
    gt_assert(cntlist != NULL);
    if (alloc_cntlist)
    {
      GT_INITBITTAB(*cntlist, *nofreads);
      n = fread(*cntlist, sizeof (GtBitsequence),
          GT_NUMOFINTSFORBITS(*nofreads), infp);
      if (n != GT_NUMOFINTSFORBITS(*nofreads))
      {
        gt_error_set(err, "contained reads file: unrecognized format");
        had_err = -1;
      }
    }
    else
    {
      /* combine using OR with existing data */
      size_t i;
      for (i = 0; i < GT_NUMOFINTSFORBITS(*nofreads); i++)
      {
        GtBitsequence value;
        n = fread(&value, sizeof (GtBitsequence), (size_t)1, infp);
        if (n != (size_t)1)
        {
          gt_error_set(err, "contained reads file: unrecognized format");
          had_err = -1;
          break;
        }
        *cntlist[i] |= value;
      }
    }
  }
  return had_err;
}

static int gt_cntlist_parse_bin(FILE *infp, bool alloc_cntlist,
    GtBitsequence **cntlist, unsigned long *nofreads, GtError *err)
{
  int had_err = gt_cntlist_parse_bin_or_bit_header(infp, nofreads, err);
  if (had_err == 0)
  {
    size_t n;
    unsigned long seqnum;
    gt_assert(cntlist != NULL);
    if (alloc_cntlist)
      GT_INITBITTAB(*cntlist, *nofreads);
    while (true)
    {
      n = fread(&seqnum, sizeof (unsigned long), (size_t)1, infp);
      if (n != (size_t)1)
      {
        if (!feof(infp))
        {
          gt_error_set(err, "contained reads file: unrecognized format");
          had_err = -1;
        }
        break;
      }
      GT_SETIBIT(*cntlist, seqnum);
    }
  }
  return had_err;
}

static int gt_cntlist_parse_ascii(FILE *infp, bool alloc_cntlist,
    GtBitsequence **cntlist, unsigned long *nofreads, GtError *err)
{
  int n;
  unsigned long seqnum;

  gt_assert(infp != NULL && nofreads != NULL && cntlist != NULL);
  /*@i1@*/ gt_error_check(err);
  n = fscanf(infp, "[n: %lu]\n", nofreads);
  if (n!=1 || *nofreads == 0)
  {
    gt_error_set(err, "contained reads file: unrecognized format");
    return -1;
  }
  if (alloc_cntlist)
    GT_INITBITTAB(*cntlist, *nofreads);
  while (true)
  {
    n = fscanf(infp, "%lu\n", &seqnum);
    if (n == EOF)
      break;
    else if (n != 1)
    {
      gt_error_set(err, "contained reads file: unrecognized format");
      return -1;
    }
    GT_SETIBIT(*cntlist, seqnum);
  }
  return 0;
}

int gt_cntlist_parse(const char *filename, bool alloc_cntlist,
    GtBitsequence **cntlist, unsigned long *nofreads, GtError *err)
{
  int c, retval = 0;
  FILE *infp;

  gt_log_log("parse contained reads list file: %s", filename);
  infp = gt_fa_fopen(filename, "rb", err);

  if (infp == NULL)
    return -1;

  c = gt_xfgetc(infp);
  switch (c)
  {
    case EOF:
      gt_error_set(err, "%s: unexpected end of file", filename);
      retval = 1;
      break;
    case GT_CNTLIST_BIN_HEADER:
      gt_log_log("contained reads list format: BIN");
      retval = gt_cntlist_parse_bin(infp, alloc_cntlist, cntlist, nofreads,
          err);
      break;
    case GT_CNTLIST_BIT_HEADER:
      gt_log_log("contained reads list format: BIT");
      retval = gt_cntlist_parse_bit(infp, alloc_cntlist, cntlist, nofreads,
          err);
      break;
    case GT_CNTLIST_ASCII_HEADER:
      gt_xungetc(c, infp);
      gt_log_log("contained reads list format: ASCII");
      retval = gt_cntlist_parse_ascii(infp, alloc_cntlist, cntlist, nofreads,
          err);
      break;
    default:
      gt_error_set(err, "%s: unrecognized format", filename);
      retval = 1;
      break;
  }
  gt_fa_fclose(infp);

  return retval;
}

unsigned long gt_cntlist_count(const GtBitsequence *cntlist,
    unsigned long nofreads)
{
  unsigned long i, counter = 0;

  for (i = 0; i < nofreads; i++)
    if ((bool)GT_ISIBITSET(cntlist, i))
      counter++;
  return counter;
}

unsigned long gt_cntlist_xload(const char *filename, GtBitsequence **cntlist,
    unsigned long expected_nofreads)
{
  int retval;
  unsigned long found_nofreads;
  GtError *err;

  if (!gt_file_exists(filename))
  {
    fprintf(stderr, "FATAL: error by loading contained reads list: "
        "file %s does not exist\n", filename);
    exit(EXIT_FAILURE);
  }

  err = gt_error_new();
  retval = gt_cntlist_parse(filename, true, cntlist, &found_nofreads, err);
  if (retval != 0)
  {
    fprintf(stderr, "FATAL: error by parsing contained reads list: %s\n",
        gt_error_get(err));
    exit(EXIT_FAILURE);
  }
  gt_error_delete(err);

  if (found_nofreads != expected_nofreads)
  {
    fprintf(stderr, "FATAL: error by parsing contained reads list: "
        "file specifies a wrong number of reads\nexpected %lu, found %lu\n",
        expected_nofreads, found_nofreads);
    exit(EXIT_FAILURE);
  }

  return gt_cntlist_count(*cntlist, found_nofreads);
}
