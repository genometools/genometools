/*
  Copyright (c) 2013 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

#include <string.h>

#include "core/fasta_api.h"
#include "core/fasta_separator.h"
#include "core/log_api.h"
#include "core/logger.h"
#include "core/ma.h"
#include "core/output_file_api.h"
#include "core/showtime.h"
#include "core/str_api.h"
#include "core/types_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "extended/condenseq.h"
#include "tools/gt_condenseq_extract.h"

typedef struct {
  GtFile           *outfp;
  GtOutputFileInfo *ofi;
  GtStr            *mode,
                   *sepchar;
  GtOption         *sepchar_opt;
  GtRange           range,
                    seqrange;
  GtUword           seq,
                    width;
  bool              verbose;
} GtCondenserExtractArguments;

static void* gt_condenseq_extract_arguments_new(void)
{
  GtCondenserExtractArguments *arguments = gt_calloc((size_t) 1,
                                                     sizeof *arguments);
  arguments->range.start =
    arguments->range.end = GT_UNDEF_UWORD;
  arguments->seqrange.start =
    arguments->seqrange.end = GT_UNDEF_UWORD;
  arguments->ofi = gt_output_file_info_new();
  arguments->mode = gt_str_new();
  arguments->sepchar = gt_str_new();
  arguments->sepchar_opt = NULL;
  return arguments;
}

static void gt_condenseq_extract_arguments_delete(void *tool_arguments)
{
  GtCondenserExtractArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_file_delete(arguments->outfp);
    gt_output_file_info_delete(arguments->ofi);
    gt_str_delete(arguments->mode);
    gt_str_delete(arguments->sepchar);
    gt_option_delete(arguments->sepchar_opt);
    gt_free(arguments);
  }
}

static GtOptionParser*
gt_condenseq_extract_option_parser_new(void *tool_arguments)
{
  GtCondenserExtractArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option,
           *optionrange,
           *optionseq,
           *optionseqrange,
           *optionmode;
  static const char *modes[] = {"fasta", "concat", NULL};
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[options] archive",
                            "Decompresses condenseq archive.");

  /* -seq */
  optionseq = gt_option_new_uword("seq",
                                  "only extract sequence identified by its "
                                  "number",
                                  &arguments->seq, GT_UNDEF_UWORD);
  gt_option_parser_add_option(op, optionseq);

  /* -seqrange */
  optionseqrange = gt_option_new_range("seqrange",
                                       "only extract (inclusive) range of "
                                       "consecutive sequences identified by "
                                       "their zero based index numbers",
                                       &arguments->seqrange, NULL);
  gt_option_parser_add_option(op, optionseqrange);
  gt_option_exclude(optionseq, optionseqrange);

  /* -range */
  optionrange = gt_option_new_range("range",
                                    "only extract (inclusive) range of zero "
                                    "based positions "
                                    "(implies option -output concat)",
                                    &arguments->range, NULL);
  gt_option_parser_add_option(op, optionrange);
  gt_option_exclude(optionseq, optionrange);
  gt_option_exclude(optionseqrange, optionrange);

  /* -output */
  optionmode = gt_option_new_choice("output",
                                    "specify output format "
                                    "(choose from fasta|concat)",
                                    arguments->mode,
                                    modes[0],
                                    modes);
  gt_option_parser_add_option(op, optionmode);
  gt_option_imply(optionrange, optionmode);

  /* -sepchar */
  option = gt_option_new_string("sepchar",
                                "specify character to print as SEPARATOR "
                                "(implies option -output concat",
                                arguments->sepchar, "|");
  gt_option_parser_add_option(op, option);
  gt_option_imply(option, optionmode);
  arguments->sepchar_opt = gt_option_ref(option);

  gt_output_file_info_register_options(arguments->ofi, op, &arguments->outfp);

  /* -width */
  option = gt_option_new_width(&arguments->width);
  gt_option_parser_add_option(op, option);

  /* -verbose */
  option = gt_option_new_bool("verbose", "Print out verbose output to stderr.",
                              &arguments->verbose, false);
  gt_option_parser_add_option(op, option);
  return op;
}

static int gt_condenseq_extract_arguments_check(int rest_argc,
                                                void *tool_arguments,
                                                GtError *err)
{
  GtCondenserExtractArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  if (rest_argc != 1) {
    had_err = -1;
    gt_error_set(err, "archive parameter is mandatory, use -help for usage");
  }
  if (!had_err &&
      arguments->range.start != GT_UNDEF_UWORD) {
    if (arguments->range.start > arguments->range.end) {
      had_err = -1;
      gt_error_set(err, "give range of positions in the form 'a b' with a "
                   "<= b");
    }
    if (!had_err &&
        strcmp(gt_str_get(arguments->mode), "fasta") == 0) {
      had_err = -1;
      gt_error_set(err, "-range is incompatible to fasta output");
    }
  }
  if (!had_err &&
      arguments->seqrange.start != GT_UNDEF_UWORD) {
    if (arguments->seqrange.start > arguments->seqrange.end) {
      had_err = -1;
      gt_error_set(err, "give range of positions in the form 'a b' with a "
                   "<= b");
    }
  }

  if (!had_err &&
      gt_option_is_set(arguments->sepchar_opt)) {
    if (strcmp(gt_str_get(arguments->mode), "fasta") == 0) {
      had_err = -1;
      gt_error_set(err, "-sepchar is incompatible with -output fasta");
    }
    if (!had_err &&
        gt_str_length(arguments->sepchar) != (GtUword) 1) {
      had_err = -1;
      gt_error_set(err, "sepchar (\"%s\" given) has to be a single character",
                   gt_str_get(arguments->sepchar));
    }
  }

  return had_err;
}

static int gt_condenseq_extract_runner(GT_UNUSED int argc,
                                       const char **argv,
                                       int parsed_args,
                                       void *tool_arguments,
                                       GtError *err)
{
  int had_err = 0;
  GtCondenserExtractArguments *arguments = tool_arguments;
  GtCondenseq *condenseq = NULL;
  GtLogger *logger = NULL;
  GtTimer *timer = NULL;

  gt_error_check(err);
  gt_assert(arguments);

  if (gt_showtime_enabled()) {
    timer = gt_timer_new_with_progress_description("load condenseq");
    gt_timer_start(timer);
  }

  logger = gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX, stderr);

  if (!had_err) {
    condenseq = gt_condenseq_new_from_file(argv[parsed_args], logger, err);
    if (condenseq == NULL) {
      had_err = -1;
    }
  }

  if (!had_err) {
    const char *buffer = NULL;
    const char *desc = NULL;
    GtUword desclen,
            seqlen,
            rend = gt_condenseq_total_length(condenseq),
            send = gt_condenseq_num_of_sequences(condenseq);
    bool concat = strcmp(gt_str_get(arguments->mode), "concat") == 0;
    /* single sequence to extract = range of length 1 */
    if (arguments->seq != GT_UNDEF_UWORD) {
      arguments->seqrange.start = arguments->seqrange.end = arguments->seq;
    }
    /* no range given at all: extract all seqs */
    if (arguments->range.start == GT_UNDEF_UWORD &&
        arguments->seqrange.start == GT_UNDEF_UWORD) {
      arguments->seqrange.start = 0;
      arguments->seqrange.end = send - 1;
    }
    /* if seqs are specified, and concat is given, switch to posrange */
    if (concat && arguments->seqrange.start != GT_UNDEF_UWORD) {
      if (arguments->seqrange.end >= send) {
        had_err = -1;
        gt_error_set(err, "range end " GT_WU " excedes number of sequences "
                     GT_WU " (ranges are zero based sequence ids)",
                     arguments->seqrange.end, send);
      }
      else {
        arguments->range.start =
          gt_condenseq_seqstartpos(condenseq, arguments->seqrange.start);
        arguments->range.end =
          gt_condenseq_seqstartpos(condenseq, arguments->seqrange.end) +
          gt_condenseq_seqlength(condenseq, arguments->seqrange.end) - 1;
      }
    }
    /* extract sequence region */
    if (!had_err && arguments->range.start != GT_UNDEF_UWORD) {
      const GtUword maxbuffsize = ((GtUword) 1) << 17; /* ~ 100000byte */
      GtUword clen,
              rstart,
              current_length = 0, i;
      const char sepchar = gt_str_get(arguments->sepchar)[0];

      if (timer)
        gt_timer_show_progress(timer, "extract sequence region(s)", stderr);

      if (arguments->range.end >= rend) {
        had_err = -1;
        gt_error_set(err, "range end " GT_WU " excedes length of sequence "
                     GT_WU " (ranges are zero based positions)",
                     arguments->range.end, rend);
      }
      if (!had_err) {
        rstart = arguments->range.start;
        rend = arguments->range.end;
        /* nextlength = gt_condenseq_seqlength(condenseq, seqnum); */
        /* seqstart = gt_condenseq_seqstartpos(condenseq, seqnum); */
        /* gt_assert(rstart >= seqstart); */
        /* nextlength -= rstart - seqstart; [> handle first seq <] */
        while (rstart <= rend) {
          GtRange cur_range;
          if (rend - rstart > maxbuffsize) {
            GtUword seqnum = gt_condenseq_pos2seqnum(condenseq,
                                                     rstart + maxbuffsize),
                    closest_sep = gt_condenseq_seqstartpos(condenseq,
                                                           seqnum) - 1;
            gt_assert(closest_sep > rstart);
            clen = closest_sep - rstart + 1;
          }
          else
            clen = rend - rstart + 1;

          cur_range.start = rstart;
          cur_range.end = rstart + clen - 1;
          buffer = gt_condenseq_extract_decoded_range(condenseq, cur_range,
                                                      sepchar);
          gt_assert(buffer != NULL);
          for (i = 0; i < clen; i++, current_length++) {
            if (arguments->width && current_length == arguments->width) {
              gt_file_xfputc('\n', arguments->outfp);
              current_length = 0;
            }
            gt_file_xfputc(buffer[i], arguments->outfp);
          }
          rstart += clen;
        }
        gt_file_xfputc('\n', arguments->outfp);
      }
    }
    else if (!had_err) { /* extract seqwise and always fasta */
      GtUword seqnum,
              sstart = arguments->seqrange.start;

      if (timer)
        gt_timer_show_progress(timer, "extract sequence(s)", stderr);
      if (arguments->seqrange.end >= send) {
        had_err = -1;
        gt_error_set(err, "range end " GT_WU " excedes number of sequences "
                     GT_WU " (ranges are zero based sequence ids)",
                     arguments->seqrange.end, send);
      }
      send = arguments->seqrange.end;
      for (seqnum = sstart;
           !had_err && seqnum <= send;
           ++seqnum) {
        buffer = gt_condenseq_extract_decoded(condenseq, &seqlen, seqnum);
        desc = gt_condenseq_description(condenseq, &desclen, seqnum);
        gt_fasta_show_entry_nt(desc, desclen,
                               buffer, seqlen,
                               arguments->width,
                               arguments->outfp);
      }
    }
  }
  if (timer)
    gt_timer_show_progress(timer, "cleanup", stderr);
  gt_condenseq_delete(condenseq);
  gt_logger_delete(logger);
  if (timer)
    gt_timer_show_progress_final(timer, stderr);
  gt_timer_delete(timer);
  return had_err;
}

GtTool* gt_condenseq_extract(void)
{
  return gt_tool_new(gt_condenseq_extract_arguments_new,
                     gt_condenseq_extract_arguments_delete,
                     gt_condenseq_extract_option_parser_new,
                     gt_condenseq_extract_arguments_check,
                     gt_condenseq_extract_runner);
}
