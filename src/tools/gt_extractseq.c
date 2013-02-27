/*
  Copyright (c) 2007-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2009 Stefan Kurtz <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/bioseq_iterator.h"
#include "core/fasta.h"
#include "core/grep_api.h"
#include "core/ma.h"
#include "core/option_api.h"
#include "core/output_file_api.h"
#include "core/unused_api.h"
#include "match/giextract.h"
#include "extended/gtdatahelp.h"
#include "tools/gt_extractseq.h"

#define FROMPOS_OPTION_STR  "frompos"
#define TOPOS_OPTION_STR    "topos"

typedef struct {
  GtStr *pattern,
        *fastakeyfile;
  unsigned long frompos,
                topos,
                width;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} ExtractSeqArguments;

static void* gt_extractseq_arguments_new(void)
{
  ExtractSeqArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->pattern = gt_str_new();
  arguments->fastakeyfile = gt_str_new();
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_extractseq_arguments_delete(void *tool_arguments)
{
  ExtractSeqArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_str_delete(arguments->fastakeyfile);
  gt_str_delete(arguments->pattern);
  gt_free(arguments);
}

static GtOptionParser* gt_extractseq_option_parser_new(void *tool_arguments)
{
  ExtractSeqArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *frompos_option, *topos_option, *match_option, *width_option,
         *fastakeyfile_option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [sequence_file(s)] | fastaindex",
                            "Extract sequences from given sequence file(s) or "
                            "fastaindex.");

  /* -frompos */
  frompos_option = gt_option_new_ulong_min(FROMPOS_OPTION_STR,
                                        "extract sequence from this position\n"
                                        "counting from 1 on",
                                        &arguments->frompos, 0, 1UL);
  gt_option_parser_add_option(op, frompos_option);

  /* -topos */
  topos_option = gt_option_new_ulong_min(TOPOS_OPTION_STR,
                                      "extract sequence up to this position\n"
                                      "counting from 1 on",
                                      &arguments->topos, 0, 1UL);
  gt_option_parser_add_option(op, topos_option);

  /* -match */
  match_option = gt_option_new_string("match", "extract all sequences whose "
                                   "description matches the given pattern.\n"
                                   "The given pattern must be a valid extended "
                                   "regular expression.", arguments->pattern,
                                   NULL);
  gt_option_parser_add_option(op, match_option);

  /* -keys */
  fastakeyfile_option = gt_option_new_filename("keys",
                                               "extract substrings for keys "
                                               "in specified file",
                                     arguments->fastakeyfile);
  gt_option_parser_add_option(op, fastakeyfile_option);

  /* -width */
  width_option = gt_option_new_width(&arguments->width);
  gt_option_parser_add_option(op, width_option);

  /* output file options */
  gt_output_file_register_options(op, &arguments->outfp, arguments->ofi);

  /* option implications */
  gt_option_imply(frompos_option, topos_option);
  gt_option_imply(topos_option, frompos_option);

  /* option exclusions */
  gt_option_exclude(frompos_option, match_option);
  gt_option_exclude(topos_option, match_option);
  gt_option_exclude(frompos_option, fastakeyfile_option);
  gt_option_exclude(match_option, fastakeyfile_option);

  gt_option_parser_set_comment_func(op, gt_gtdata_show_help, NULL);
  return op;
}

static int gt_extractseq_arguments_check(GT_UNUSED int argc,
                                         void *tool_arguments, GtError *err)
{
  ExtractSeqArguments *arguments = tool_arguments;
  gt_error_check(err);
  gt_assert(arguments);
  if (arguments->frompos > arguments->topos) {
    gt_error_set(err,
              "argument to option '-%s' must be <= argument to option '-%s'",
              FROMPOS_OPTION_STR, TOPOS_OPTION_STR);
    return -1;
  }
  return 0;
}

static int extractseq_pos(GtFile *outfp, GtBioseq *bs,
                          unsigned long frompos, unsigned long topos,
                          unsigned long width, GtError *err)
{
  int had_err = 0;
  GtStr *buf;
  char *out = NULL;
  unsigned long accupos = 0,
                newstartpos = 0,
                len = topos - frompos + 1,
                i = 0;
  gt_assert(bs);
  gt_error_check(err);

  if (frompos > gt_bioseq_get_total_length(bs)
        || topos > gt_bioseq_get_total_length(bs)) {
    gt_error_set(err, "invalid position pair %lu-%lu one value is larger than "
                      "sequence length %lu", frompos, topos,
                      gt_bioseq_get_total_length(bs));
    return -1;
  }
  frompos--; topos--;
  buf = gt_str_new();

  /* look for beginning of sequence */
  while (accupos + gt_bioseq_get_sequence_length(bs, i) <= frompos
            && i < gt_bioseq_number_of_sequences(bs)) {
    accupos += gt_bioseq_get_sequence_length(bs, i);
    i++;
  }
  if (i == 0) {
    newstartpos = frompos;
    accupos = frompos;
  } else {
    gt_assert(accupos > 0);
    newstartpos = frompos - accupos;
  }

  /* do we need to cross a sequence boundary to print the full output? */
  if (len <= gt_bioseq_get_sequence_length(bs, i) - newstartpos) {
    /* no, just print */
    out = gt_bioseq_get_sequence_range(bs, i, newstartpos,
                                       newstartpos + len - 1);
    gt_str_append_cstr_nt(buf, out, len);
    gt_free(out);
  } else {
    /* yes, first output the part on this sequence... */
    unsigned long restlen = gt_bioseq_get_sequence_length(bs, i) - newstartpos,
                  restfulllen = topos - accupos + 1;
    out = gt_bioseq_get_sequence_range(bs, i, newstartpos,
                                       newstartpos + restlen - 1);
    restfulllen -= restlen;
    gt_str_append_cstr_nt(buf, out, restlen);
    gt_free(out);

    i++;
    /* ...then determine whether we need to output full seqs in between... */
    while (restfulllen > gt_bioseq_get_sequence_length(bs, i)
             && i < gt_bioseq_number_of_sequences(bs) - 1) {
      unsigned long thislen = gt_bioseq_get_sequence_length(bs, i);
      out = gt_bioseq_get_sequence_range(bs, i, 0, thislen - 1);
      gt_str_append_cstr_nt(buf, out, thislen);
      gt_free(out);
      restfulllen -= thislen;
      i++;
    }

    /* ...then output the last sequence */
    if (restfulllen > 0) {
      out = gt_bioseq_get_sequence_range(bs, i, 0, restfulllen - 1);
      gt_str_append_cstr_nt(buf, out, restfulllen);
      gt_free(out);
    }
  }

  gt_fasta_show_entry(NULL, gt_str_get(buf), gt_str_length(buf), width, outfp);
  gt_str_delete(buf);
  return had_err;
}

static int extractseq_match(GtFile *outfp, GtBioseq *bs,
                            const char *pattern, unsigned long width,
                            GtError *err)
{
  const char *desc;
  unsigned long i;
  bool match;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(bs && pattern);

  for (i = 0; !had_err && i < gt_bioseq_number_of_sequences(bs); i++) {
    desc = gt_bioseq_get_description(bs, i);
    gt_assert(desc);
    had_err = gt_grep(&match, pattern, desc, err);
    if (!had_err && match) {
      char *out = gt_bioseq_get_sequence(bs, i);
      gt_fasta_show_entry(desc, out, gt_bioseq_get_sequence_length(bs, i),
                          width, outfp);
      gt_free(out);
    }
  }

  return had_err;
}

static int process_fastakeyfile(GtStr *fastakeyfile, int argc,
                                const char **argv, unsigned long width,
                                GtFile *outfp, GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  gt_assert(gt_str_length(fastakeyfile));

  if (argc == 0) {
    gt_error_set(err,"option -keys requires at least one file argument");
    had_err = -1;
  }

  if (!had_err)
  {
    const char *indexname = argv[0];

    if (argc == 1 && gt_deskeysfileexists(indexname))
    {
      if (gt_extractkeysfromfastaindex(indexname,fastakeyfile,width,err) != 0)
      {
        had_err = -1;
      }
    } else
    {
      GtStrArray *referencefiletab;
      int i;

      referencefiletab = gt_str_array_new();
      for (i = 0; i < argc; i++)
      {
        gt_str_array_add_cstr(referencefiletab, argv[i]);
      }
      if (gt_extractkeysfromfastafile(true, outfp, width, fastakeyfile,
                                      referencefiletab, err) != 1)
      {
        had_err = -1;
      }
      gt_str_array_delete(referencefiletab);
    }
  }
  return had_err;
}

static int gt_extractseq_runner(int argc, const char **argv, int parsed_args,
                                void *tool_arguments, GtError *err)
{
  ExtractSeqArguments *arguments = tool_arguments;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);
  if (gt_str_length(arguments->fastakeyfile)) {
    had_err = process_fastakeyfile(arguments->fastakeyfile, argc - parsed_args,
                            argv + parsed_args, arguments->width,
                            arguments->outfp, err);
  }
  else {
    GtBioseqIterator *bsi;
    GtBioseq *bs;
    bsi = gt_bioseq_iterator_new(argc - parsed_args, argv + parsed_args);
    while (!had_err &&
           !(had_err = gt_bioseq_iterator_next(bsi, &bs, err)) && bs) {
      if (arguments->frompos) {
        had_err = extractseq_pos(arguments->outfp, bs, arguments->frompos,
                                 arguments->topos, arguments->width, err);
      }
      else {
        had_err = extractseq_match(arguments->outfp, bs,
                                   gt_str_get(arguments->pattern),
                                   arguments->width, err);
      }
      gt_bioseq_delete(bs);
    }
    gt_bioseq_iterator_delete(bsi);
  }
  return had_err;
}

GtTool* gt_extractseq(void)
{
  return gt_tool_new(gt_extractseq_arguments_new,
                     gt_extractseq_arguments_delete,
                     gt_extractseq_option_parser_new,
                     gt_extractseq_arguments_check,
                     gt_extractseq_runner);
}
