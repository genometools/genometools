/*
  Copyright (c) 2009-2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2009-2011 Center for Bioinformatics, University of Hamburg

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

#include <errno.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include "core/str.h"
#include "core/assert_api.h"
#include "core/file.h"
#include "core/log_api.h"
#include "core/parseutils.h"
#include "core/splitter.h"
#include "core/ensure.h"
#include "core/minmax.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "match/rdj-spmlist.h"
/* unit test: */
#include "match/rdj-spmproc.h"
#include "match/rdj-ensure-output.h"

/* ---------------- Readjoiner bin format ---------------- */

#ifndef NDEBUG
#define GT_SPMLIST_ASSERT_CAST_SAFE(N, FROM, TO, TOMAX)\
  if (sizeof (FROM) > sizeof (TO)) \
    gt_assert(N < (FROM)TOMAX)
#else
#define GT_SPMLIST_ASSERT_CAST_SAFE(N, FROM, TO, TOMAX)
#endif

/*
 * INT1: suffix_seqnum
 * INT2: prefix_seqnum
 * INT3: length << 2 + suffixseq_direct << 1 + prefixseq_direct
 *
 * */

#if defined(_LP64) || (BITS == 32)
#define GT_RDJ_LENGTHASSERTION(BITS)\
        gt_assert(length <= (UINT ## BITS ## _MAX >> 2))
#else
#define GT_RDJ_LENGTHASSERTION(BITS) /* Nothing */
#endif

#define DEFINE_GT_SPMLIST_BIN_FORMAT(BITS)\
void gt_spmlist_write_header_bin ## BITS(FILE *file)\
{\
  gt_xfputc((int)GT_SPMLIST_BIN ## BITS, file);\
}\
void gt_spmproc_show_bin ## BITS(unsigned long suffix_seqnum,\
    unsigned long prefix_seqnum, unsigned long length, bool suffixseq_direct,\
    bool prefixseq_direct, void *file)\
{\
  uint ## BITS ## _t spmdata[3];\
  length <<= 2;\
  GT_RDJ_LENGTHASSERTION(BITS);\
  if (suffixseq_direct)\
    length |= 2;\
  if (prefixseq_direct)\
    length |= 1;\
  GT_SPMLIST_ASSERT_CAST_SAFE(suffix_seqnum, unsigned long, uint ## BITS ## _t,\
      UINT ## BITS ## _MAX);\
  spmdata[0] = (uint ## BITS ## _t) suffix_seqnum;\
  GT_SPMLIST_ASSERT_CAST_SAFE(prefix_seqnum, unsigned long, uint ## BITS ## _t,\
      UINT ## BITS ## _MAX);\
  spmdata[1] = (uint ## BITS ## _t) prefix_seqnum;\
  spmdata[2] = (uint ## BITS ## _t) length;\
  /*@ignore@*/\
  gt_xfwrite(&spmdata, sizeof (uint ## BITS ## _t), (size_t)3, (FILE*)file);\
  /*@end@*/\
}\
static int gt_spmlist_parse_bin ## BITS(FILE *file, unsigned long min_length,\
    GtSpmproc processoverlap, void *data, GtError *err)\
{\
  int had_err = 0;\
  size_t retval;\
  uint ## BITS ## _t spmdata[3];\
  unsigned long length;\
  bool suffixseq_direct, prefixseq_direct;\
  while (!feof(file))\
  {\
    retval = fread(&spmdata, sizeof (uint ## BITS ## _t), (size_t)3, file);\
    if (retval == 0 && feof(file))\
      break;\
    if (retval != (size_t)3)\
    {\
      had_err = -1;\
      gt_log_log("retval: %lu", (unsigned long)retval);\
      gt_error_set(err, "SPM binary file error: %s", feof(file) ?\
          "premature EOF" : strerror(errno));\
    }\
    GT_SPMLIST_ASSERT_CAST_SAFE(spmdata[2] >> 2, uint ## BITS ## _t,\
        unsigned long, ULONG_MAX);\
    length = (unsigned long)(spmdata[2] >> 2);\
    suffixseq_direct = (spmdata[2] & 2) != 0;\
    prefixseq_direct = (spmdata[2] & 1) != 0;\
    GT_SPMLIST_ASSERT_CAST_SAFE(spmdata[0], uint ## BITS ## _t,\
        unsigned long, ULONG_MAX);\
    GT_SPMLIST_ASSERT_CAST_SAFE(spmdata[1], uint ## BITS ## _t,\
        unsigned long, ULONG_MAX);\
    if (had_err == 0 && length >= min_length)\
      processoverlap((unsigned long)spmdata[0], (unsigned long)spmdata[1],\
          length, suffixseq_direct, prefixseq_direct, data);\
  }\
  return had_err;\
}

DEFINE_GT_SPMLIST_BIN_FORMAT(32);
DEFINE_GT_SPMLIST_BIN_FORMAT(64);

/* ---------------- Plain text format ---------------- */

/*@notfunction@*/
#define GT_SPMLIST_EOFERROR\
  do {\
    gt_error_set(err, "unexpected end of file");\
    return -1;\
  } while (false)

static inline int parse_plusminus(bool *destination, const char *source)
{
  if (strlen(source) != (size_t)1)
    return -1;
  switch (source[0])
  {
    case '+':
      *destination = true;
      return 0;
    case '-':
      *destination = false;
      return 0;
    default:
      return -1;
  }
}

#define GT_SPMLIST_PARSE(NR, F, PTR)\
  if (!had_err) \
  {\
    gt_assert(tokens != NULL);\
    had_err = F((PTR), tokens[NR]);\
    if (had_err) \
      gt_error_set(err, "Token %i unrecognized", NR);\
  }

static inline int parse_line(GtStr *s, unsigned long min_length,
    GtSpmproc proc_e, GtSpmprocA proc_a, void *data, GtError *err)
{
  int had_err = 0;
  GtSplitter *splitter;
  char **tokens = NULL;
  unsigned long suffix_seqnum = 0, prefix_seqnum = 0, suffix_length = 0,
                prefix_length, unit_edist;
  bool suffixseq_direct = true, prefixseq_direct = true, exact;

  exact = (proc_e != NULL) ? true : false;
  splitter = gt_splitter_new();
  gt_splitter_split(splitter, gt_str_get(s), gt_str_length(s), ' ');
  if (gt_splitter_size(splitter) != (exact ? 5UL : 7UL))
  {
    gt_error_set(err, "Wrong number of tokens");
    had_err = -1;
  }
  if (!had_err) tokens = gt_splitter_get_tokens(splitter);
  GT_SPMLIST_PARSE(0, gt_parse_ulong, &suffix_seqnum);
  GT_SPMLIST_PARSE(1, parse_plusminus, &suffixseq_direct);
  GT_SPMLIST_PARSE(2, gt_parse_ulong, &prefix_seqnum);
  GT_SPMLIST_PARSE(3, parse_plusminus, &prefixseq_direct);
  GT_SPMLIST_PARSE(4, gt_parse_ulong, &suffix_length);
  if (!exact)
  {
    GT_SPMLIST_PARSE(5, gt_parse_ulong, &prefix_length);
    GT_SPMLIST_PARSE(6, gt_parse_ulong, &unit_edist);
  }
  if (!had_err)
  {
    if (exact && suffix_length >= min_length)
    {
      gt_assert(proc_e != NULL);
      proc_e(suffix_seqnum, prefix_seqnum, suffix_length,
          suffixseq_direct, prefixseq_direct, data);
    }
    else if (!exact && (suffix_length >= min_length ||
          prefix_length >= min_length))
    {
      gt_assert(proc_a != NULL);
      proc_a(suffix_seqnum, prefix_seqnum, suffix_length, prefix_length,
          unit_edist, suffixseq_direct, prefixseq_direct, data);
    }
  }
  gt_splitter_delete(splitter);
  return had_err;
}

static int gt_spmlist_parse_ascii_generic(GtFile *infp,
    unsigned long min_length, GtSpmproc proc_e, GtSpmprocA proc_a, void *data,
    GtError *err)
{
  int had_err = 0;
  GtStr *line;
  gt_error_check(err);
  line = gt_str_new();
  while (!had_err && gt_str_read_next_line_generic(line, infp) != EOF)
  {
    had_err = parse_line(line, min_length, proc_e, proc_a, data, err);
    gt_str_reset(line);
  }
  gt_str_delete(line);
  return had_err;
}

static inline int gt_spmlist_parse_ascii(GtFile *infp, unsigned long min_length,
    GtSpmproc processoverlap, void *data, GtError *err)
{
  return gt_spmlist_parse_ascii_generic(infp, min_length, processoverlap, NULL,
      data, err);
}

int gt_spmlist_parse_ascii_approx(const char* filename,
    unsigned long min_length, GtSpmprocA processoverlap, void *data,
    GtError *err)
{
  int retval = 0;
  GtFile *infp;

  /*@i1@*/ gt_error_check(err);
  infp = gt_file_new(filename, "r", err);
  if (infp == NULL) return -1;
  retval = gt_spmlist_parse_ascii_generic(infp, min_length, NULL,
      processoverlap, data, err);
  gt_file_delete(infp);
  return retval;
}

int gt_spmlist_parse(const char* filename, unsigned long min_length,
    GtSpmproc processoverlap, void *data, GtError *err)
{
  int c, retval = 0;
  FILE *file;
  GtFile *infp;

  file = gt_fa_fopen(filename, "rb", err);
  if (file == NULL)
    return -1;
  infp = gt_file_new_from_fileptr(file);
  if (infp == NULL)
    return -1;

  c = gt_file_xfgetc(infp);
  switch (c)
  {
    case EOF:
      gt_error_set(err, "%s: file is empty", filename);
      retval = -1;
      break;
    case GT_SPMLIST_BIN32:
      gt_log_log("Spm file %s format: readjoiner-bin32", filename);
      retval = gt_spmlist_parse_bin32(file, min_length, processoverlap, data,
          err);
      break;
    case GT_SPMLIST_BIN64:
      gt_log_log("Spm file %s format: readjoiner-bin64", filename);
      retval = gt_spmlist_parse_bin64(file, min_length, processoverlap, data,
          err);
      break;
    default:
      gt_file_unget_char(infp, c);
      gt_log_log("Spm file %s format: readjoiner-text", filename);
      retval = gt_spmlist_parse_ascii(infp, min_length, processoverlap, data,
          err);
  }
  gt_file_delete(infp);
  return retval;
}

/* -------------------------- unit tests -------------------------- */

static inline int parse_plusminus_unit_test(GtError *err)
{
  int had_err = 0;
  bool destination;
  gt_ensure(had_err, parse_plusminus(&destination, "+") == 0);
  gt_ensure(had_err, destination);
  gt_ensure(had_err, parse_plusminus(&destination, "-") == 0);
  gt_ensure(had_err, !destination);
  gt_ensure(had_err, parse_plusminus(&destination, "xy") != 0);
  gt_ensure(had_err, parse_plusminus(&destination, "x") != 0);
  gt_ensure(had_err, parse_plusminus(&destination, "") != 0);
  return had_err;
}

struct GtSpmParseExactResult
{
  unsigned long suffix_seqnum, prefix_seqnum, length;
  bool suffixseq_direct, prefixseq_direct;
};

struct GtSpmParseApproxResult
{
  unsigned long suffix_seqnum, prefix_seqnum,
                suffix_length, prefix_length, unit_edist;
  bool suffixseq_direct, prefixseq_direct;
};

static void gt_spmlist_test_save(unsigned long suffix_seqnum,
    unsigned long prefix_seqnum, unsigned long length, bool suffixseq_direct,
    bool prefixseq_direct, void* data)
{
  struct GtSpmParseExactResult *r = data;
  r->suffix_seqnum = suffix_seqnum;
  r->prefix_seqnum = prefix_seqnum;
  r->length = length;
  r->suffixseq_direct = suffixseq_direct;
  r->prefixseq_direct = prefixseq_direct;
}

static void gt_spmlist_test_save_a(unsigned long suffix_seqnum,
    unsigned long prefix_seqnum, unsigned long suffix_length,
    unsigned long prefix_length, unsigned long unit_edist,
    bool suffixseq_direct, bool prefixseq_direct, void* data)
{
  struct GtSpmParseApproxResult *r = data;
  r->suffix_seqnum = suffix_seqnum;
  r->prefix_seqnum = prefix_seqnum;
  r->suffix_length = suffix_length;
  r->prefix_length = prefix_length;
  r->unit_edist = unit_edist;
  r->suffixseq_direct = suffixseq_direct;
  r->prefixseq_direct = prefixseq_direct;
}

#define GT_SPMPARSE_TEST_LINE_ERR(LINE)\
  do {\
    line = gt_str_new_cstr(LINE);\
    gt_ensure(had_err, parse_line(line, 0, gt_spmlist_test_save, NULL, &r,\
          parse_err) != 0);\
    gt_ensure(had_err, gt_error_is_set(parse_err));\
    gt_error_unset(parse_err);\
    gt_str_delete(line);\
  } while (false)

#define GT_SPMPARSE_TEST_LINE_OK(LINE,SN,PN,L,SD,PD)\
  do {\
    line = gt_str_new_cstr(LINE);\
    gt_ensure(had_err, parse_line(line, 0, gt_spmlist_test_save, NULL, &r,\
          parse_err) == 0);\
    gt_ensure(had_err, !gt_error_is_set(parse_err));\
    gt_ensure(had_err, r.suffix_seqnum == (SN));\
    gt_ensure(had_err, r.prefix_seqnum == (PN));\
    gt_ensure(had_err, r.length == (L));\
    gt_ensure(had_err, r.suffixseq_direct == (SD));\
    gt_ensure(had_err, r.prefixseq_direct == (PD));\
    gt_str_delete(line);\
  } while (false)

#define GT_SPMPARSE_TEST_LINE_A_ERR(LINE)\
  do {\
    line = gt_str_new_cstr(LINE);\
    gt_ensure(had_err, parse_line(line, 0, NULL, gt_spmlist_test_save_a, &r_a,\
          parse_err) != 0);\
    gt_ensure(had_err, gt_error_is_set(parse_err));\
    gt_error_unset(parse_err);\
    gt_str_delete(line);\
  } while (false)

#define GT_SPMPARSE_TEST_LINE_A_OK(LINE,SN,PN,SL,PL,UE,SD,PD)\
  do {\
    line = gt_str_new_cstr(LINE);\
    gt_ensure(had_err, parse_line(line, 0, NULL, gt_spmlist_test_save_a, &r_a,\
          parse_err) == 0);\
    gt_ensure(had_err, !gt_error_is_set(parse_err));\
    gt_ensure(had_err, r_a.suffix_seqnum == (SN));\
    gt_ensure(had_err, r_a.prefix_seqnum == (PN));\
    gt_ensure(had_err, r_a.suffix_length == (SL));\
    gt_ensure(had_err, r_a.prefix_length == (PL));\
    gt_ensure(had_err, r_a.unit_edist == (UE));\
    gt_ensure(had_err, r_a.suffixseq_direct == (SD));\
    gt_ensure(had_err, r_a.prefixseq_direct == (PD));\
    gt_str_delete(line);\
  } while (false)

static inline int parse_line_unit_test(GtError* err)
{
  int had_err = 0;
  GtStr *line;
  GtError *parse_err;
  struct GtSpmParseExactResult r;
  struct GtSpmParseApproxResult r_a;

  gt_error_check(err);
  parse_err = gt_error_new();
  GT_SPMPARSE_TEST_LINE_ERR("1 2 + 3 4");
  GT_SPMPARSE_TEST_LINE_ERR("4 5 1");
  GT_SPMPARSE_TEST_LINE_ERR("1 1 1 1");
  GT_SPMPARSE_TEST_LINE_ERR("1x 2 + 1");
  GT_SPMPARSE_TEST_LINE_OK("1 + 2 - 3", 1UL, 2UL, 3UL, true, false);
  GT_SPMPARSE_TEST_LINE_OK("2 - 1 + 3", 2UL, 1UL, 3UL, false, true);
  GT_SPMPARSE_TEST_LINE_OK("4 + 5 + 6", 4UL, 5UL, 6UL, true, true);
  GT_SPMPARSE_TEST_LINE_A_ERR("1 + 2 + 3");
  GT_SPMPARSE_TEST_LINE_A_ERR("4 + 5 + 1 5 4 5");
  GT_SPMPARSE_TEST_LINE_A_ERR("1 1 1 1 2 3");
  GT_SPMPARSE_TEST_LINE_A_ERR("1x + 2 + 1 2 3");
  GT_SPMPARSE_TEST_LINE_A_OK("1 + 2 - 3 4 1", 1UL, 2UL, 3UL, 4UL, 1UL, true,
      false);
  GT_SPMPARSE_TEST_LINE_A_OK("2 - 1 + 3 4 1", 2UL, 1UL, 3UL, 4UL, 1UL, false,
      true);
  GT_SPMPARSE_TEST_LINE_A_OK("4 + 5 + 6 7 1", 4UL, 5UL, 6UL, 7UL, 1UL, true,
      true);
  gt_error_delete(parse_err);
  return had_err;
}

static int gt_spmlist_parse_unit_test(GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  had_err = parse_plusminus_unit_test(err);
  if (!had_err) had_err = parse_line_unit_test(err);
  return had_err;
}

void gt_spmproc_show_ascii(unsigned long suffix_seqnum,
    unsigned long prefix_seqnum,
    unsigned long length,
    bool suffixseq_direct,
    bool prefixseq_direct,
    void *data)
{
  GtFile *file = data;
  gt_file_xprintf(file, "%lu %s %lu %s %lu\n", suffix_seqnum,
      suffixseq_direct ? "+" : "-", prefix_seqnum,
      prefixseq_direct ? "+" : "-", length);
}

void gt_spmproc_a_show_ascii(unsigned long suffix_seqnum,
    unsigned long prefix_seqnum,
    unsigned long suffix_length,
    unsigned long prefix_length,
    unsigned long unit_edist,
    bool suffixseq_direct,
    bool prefixseq_direct,
    void *data)
{
  GtFile *file = data;
  gt_file_xprintf(file, "%lu %s %lu %s %lu %lu %lu\n", suffix_seqnum,
      suffixseq_direct ? "+" : "-", prefix_seqnum,
      prefixseq_direct ? "+" : "-", suffix_length, prefix_length,
      unit_edist);
}

/* ---------------------- Unit Test ---------------------- */

static void spmproc_show_caller(GtSpmproc proc, void *data)
{
  proc(0UL, 1UL, 10UL, true, true, data);
  proc(0UL, 2UL, 10UL, true, false, data);
  proc(2UL, 0UL, 10UL, false, true, data);
  proc(1UL, 0UL, 10UL, false, true, data);
}

static void spmproc_a_show_caller(GtSpmprocA proc, void *data)
{
  proc(0UL, 1UL, 10UL, 11UL, 1UL, true, true, data);
  proc(0UL, 2UL, 10UL, 10UL, 0UL, true, false, data);
  proc(2UL, 0UL, 10UL, 12UL, 2UL, false, true, data);
  proc(1UL, 0UL, 11UL, 11UL, 0UL, false, true, data);
}

static int gt_spmproc_show_unit_test(GtError *err)
{
  int had_err = 0;
  GT_ENSURE_OUTPUT_DECLARE(100);

  gt_error_check(err);

  GT_ENSURE_OUTPUT(
      spmproc_show_caller(gt_spmproc_show_ascii, outfp),
      "0 + 1 + 10\n0 + 2 - 10\n2 - 0 + 10\n1 - 0 + 10\n");

  if (!had_err)
  {
    GT_ENSURE_OUTPUT(
        spmproc_a_show_caller(gt_spmproc_a_show_ascii, outfp),
        "0 + 1 + 10 11 1\n0 + 2 - 10 10 0\n"
        "2 - 0 + 10 12 2\n1 - 0 + 11 11 0\n");
  }

  return had_err;
}

int gt_spmlist_unit_test(GtError *err)
{
  int had_err = 0;
  had_err = gt_spmproc_show_unit_test(err);
  if (had_err == 0)
    had_err = gt_spmlist_parse_unit_test(err);
  return had_err;
}
