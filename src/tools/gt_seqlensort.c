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

#include "core/ma.h"
#include "core/str_array.h"
#include "core/unused_api.h"
#include "match/reads2twobit.h"
#include "tools/gt_seqlensort.h"

typedef struct {
  GtStr  *indexname;
  GtStrArray  *db;
} GtSeqlensortArguments;

static void* gt_seqlensort_arguments_new(void)
{
  GtSeqlensortArguments *arguments = gt_malloc(sizeof (*arguments));
  arguments->indexname = gt_str_new();
  arguments->db = gt_str_array_new();
  return arguments;
}

static void gt_seqlensort_arguments_delete(void *tool_arguments)
{
  GtSeqlensortArguments *arguments = tool_arguments;
  gt_str_delete(arguments->indexname);
  gt_str_array_delete(arguments->db);
  gt_free(arguments);
}

static GtOptionParser* gt_seqlensort_option_parser_new(void *tool_arguments)
{
  GtSeqlensortArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *indexname_option, *db_option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("-db <fas ...> [-indexname ...]",
      "Encode DNA MultiFasta sequences (with no wildcards) "
      "in GtEncseq format, sorting the sequences by length.");

  /* -indexname */
  indexname_option = gt_option_new_string("indexname",
      "specify the indexname to use\n"
      "default: first argument of -db option",
      arguments->indexname, NULL);
  gt_option_hide_default(indexname_option);
  gt_option_parser_add_option(op, indexname_option);

  /* -db */
  db_option = gt_option_new_filename_array("db",
      "name of input MultiFasta file(s)", arguments->db);
  gt_option_hide_default(db_option);
  gt_option_is_mandatory(db_option);
  gt_option_parser_add_option(op, db_option);

  return op;
}

#define GT_SEQLENSORT_SEQLEN(SEQNUM, SEPPOS) \
  ((SEQNUM) == 0 \
    ? ((SEPPOS)[SEQNUM]) \
    : ((SEPPOS)[SEQNUM] - ((SEPPOS)[(SEQNUM) - 1UL] + 1UL)))

static int gt_seqlensort_cmp(const void *a, const void *b, void *data)
{
  const unsigned long seqnum_a = *(const unsigned long*)a,
                      seqnum_b = *(const unsigned long*)b;
  unsigned long seqlen_a = GT_SEQLENSORT_SEQLEN(seqnum_a, (unsigned long*)data),
                seqlen_b = GT_SEQLENSORT_SEQLEN(seqnum_b, (unsigned long*)data);
  if (seqlen_a == seqlen_b)
    return (int)((seqnum_a > seqnum_b) - (seqnum_a < seqnum_b));
  return (int)((seqlen_a > seqlen_b) - (seqlen_a < seqlen_b));
}

static int gt_seqlensort_runner(GT_UNUSED int argc,
    GT_UNUSED const char **argv, GT_UNUSED int parsed_args,
    void *tool_arguments, GtError *err)
{
  GtSeqlensortArguments *arguments = tool_arguments;
  int had_err = 0;
  unsigned long i;
  GtReads2Twobit *r2t = NULL;

  gt_error_check(err);
  gt_assert(arguments);

  if (gt_str_length(arguments->indexname) == 0)
  {
    gt_str_append_cstr(arguments->indexname,
        gt_str_array_get(arguments->db, 0));
  }

  r2t = gt_reads2twobit_new(arguments->indexname);

  for (i = 0; i < gt_str_array_size(arguments->db) && !had_err; i++)
    had_err = gt_reads2twobit_add_library(r2t, gt_str_array_get_str(
          arguments->db, i), err);

  if (!had_err)
    had_err = gt_reads2twobit_encode(r2t, err);

  if (!had_err)
  {
    if (gt_reads2twobit_seqlen_eqlen(r2t) == 0)
    {
      gt_reads2twobit_sort(r2t, gt_seqlensort_cmp,
          gt_reads2twobit_export_seppos(r2t));
    }
    had_err = gt_reads2twobit_write_encseq(r2t, err);
  }

  gt_reads2twobit_delete(r2t);
  return had_err;
}

GtTool* gt_seqlensort(void)
{
  return gt_tool_new(gt_seqlensort_arguments_new,
                  gt_seqlensort_arguments_delete,
                  gt_seqlensort_option_parser_new,
                  NULL,
                  gt_seqlensort_runner);
}
