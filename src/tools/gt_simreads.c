/*
  Copyright (c) 2009-2010 Giorgio Gonnella <ggonnell@yahoo.it>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#include "core/str.h"
#include "core/disc_distri_api.h"
#include "core/encseq.h"
#include "core/fasta.h"
#include "core/fa.h"
#include "core/file_api.h"
#include "core/fileutils.h"
#include "core/logger.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/option_api.h"
#include "core/output_file_api.h"
#include "core/progressbar.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "tools/gt_simreads.h"

typedef struct {
  unsigned long num, minlen, maxlen, coverage;
  bool show_progressbar, verbose, ds, dl, singlestrand;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
  GtStr *dsfilename, *dlfilename, *distlen_filename;
  FILE *distlen_file;
} GtSimreadsArguments;

static void* gt_simreads_arguments_new(void)
{
  GtSimreadsArguments *arguments = gt_calloc((size_t)1, sizeof *arguments);
  arguments->ofi = gt_output_file_info_new();
  arguments->distlen_filename = gt_str_new();
  arguments->distlen_file = NULL;
  arguments->dsfilename = gt_str_new();
  arguments->dlfilename = gt_str_new();
  return arguments;
}

static void gt_simreads_arguments_delete(void *tool_arguments)
{
  GtSimreadsArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  if (arguments->distlen_file != NULL)
    gt_fa_fclose(arguments->distlen_file);
  gt_str_delete(arguments->distlen_filename);
  gt_str_delete(arguments->dsfilename);
  gt_str_delete(arguments->dlfilename);
  gt_free(arguments);
}

static GtOptionParser* gt_simreads_option_parser_new(void *tool_arguments)
{
  GtSimreadsArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *num_option, *coverage_option, *len_option,
           *minlen_option, *maxlen_option, *p_option, *v_option,
           *dl_option, *ds_option, *ss_option, *distlen_option;

  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] <encseq>",
                         "Simulate sequencing reads from "
                         "random positions in the input sequence(s).");

  /* -num */
  num_option = gt_option_new_ulong_min("num", "desired number of reads",
                                       &arguments->num, GT_UNDEF_ULONG, 1UL);
  gt_option_parser_add_option(op, num_option);

  /* -coverage */
  coverage_option = gt_option_new_ulong_min("coverage",
                                            "desired coverage of the reads",
                                            &arguments->coverage,
                                            GT_UNDEF_ULONG, 1UL);
  gt_option_parser_add_option(op, coverage_option);

  /* -len */
  len_option = gt_option_new_ulong_min("len", "fixed read length",
                                               &arguments->minlen,
                                               GT_UNDEF_ULONG, 1UL);
  gt_option_parser_add_option(op, len_option);

  /* -minlen */
  minlen_option = gt_option_new_ulong_min("minlen",
                                          "minimal read length",
                                          &arguments->minlen,
                                          GT_UNDEF_ULONG, 1UL);
  gt_option_parser_add_option(op, minlen_option);

  /* -maxlen */
  maxlen_option = gt_option_new_ulong_min("maxlen",
                                          "maximal read length",
                                          &arguments->maxlen,
                                          GT_UNDEF_ULONG, 1UL);
  gt_option_parser_add_option(op, maxlen_option);

  /* -distlen */
  distlen_option = gt_option_new_string("distlen",
      "use read length distribution file\n"
      "(in the output format of the seqstat tool)",
      arguments->distlen_filename, NULL);
  gt_option_parser_add_option(op, distlen_option);

  /* output file options */;
  gt_output_file_register_options(op, &arguments->outfp, arguments->ofi);

  /* -p */
  p_option = gt_option_new_bool("p", "show a progress bar",
                                &arguments->show_progressbar, false);
  gt_option_parser_add_option(op, p_option);

  /* -v */
  v_option = gt_option_new_bool("v", "be verbose",
                                &arguments->verbose, false);
  gt_option_parser_add_option(op, v_option);

  /* -ds */
  ds_option = gt_option_new_string("ds", "output distribution of starting"
                                       " positions to specified file",
                                   arguments->dsfilename, NULL);
  gt_option_parser_add_option(op, ds_option);
  gt_option_is_extended_option(ds_option);

  /* -dl */
  dl_option = gt_option_new_string("dl", "output distribution of read lengths"
                                       " to specified file",
                                   arguments->dlfilename, NULL);
  gt_option_parser_add_option(op, dl_option);
  gt_option_is_extended_option(dl_option);

  /* -ss */
  ss_option = gt_option_new_bool("ss", "simulate reads in forward direction "
                                 "only", &arguments->singlestrand, false);
  gt_option_parser_add_option(op, ss_option);

  gt_option_is_mandatory_either(num_option, coverage_option);
  gt_option_is_mandatory_either_3(len_option, minlen_option, distlen_option);
  gt_option_imply(minlen_option, maxlen_option);
  gt_option_imply(maxlen_option, minlen_option);
  gt_option_exclude(maxlen_option, len_option);
  gt_option_exclude(maxlen_option, distlen_option);
  gt_option_exclude(minlen_option, len_option);
  gt_option_exclude(minlen_option, distlen_option);
  gt_option_exclude(len_option, distlen_option);
  gt_option_exclude(len_option, dl_option);

  gt_option_parser_set_min_max_args(op, 1U, 1U);

  gt_option_parser_set_mail_address(op,"<gonnella@zbh.uni-hamburg.de>");

  return op;
}

static int gt_simreads_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments, GtError *err)
{
  GtSimreadsArguments *arguments = tool_arguments;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  arguments->dl = (gt_str_length(arguments->dlfilename) > 0);
  arguments->ds = (gt_str_length(arguments->dsfilename) > 0);

  if ((arguments->maxlen != GT_UNDEF_ULONG) &&
      (arguments->minlen > arguments->maxlen)) {
    gt_error_set(err,
      "argument to option '-minlen' must be <= argument to option '-maxlen'");
    had_err = -1;
  }

  /* cannot use imply with outputfile_register_options registered options: */
  if (!had_err && arguments->show_progressbar && !arguments->outfp)
  {
    gt_error_set(err, "option \"-p\" requires option \"-o\"");
    had_err = -1;
  }

  if (!had_err && gt_str_length(arguments->distlen_filename) > 0)
  {
    if (gt_file_exists(gt_str_get(arguments->distlen_filename)))
    {
      arguments->distlen_file =
        gt_fa_fopen(gt_str_get(arguments->distlen_filename), "rb", err);
      if (arguments->distlen_file == NULL)
        had_err = -1;
    }
    else
    {
      gt_error_set(err, "file %s not found",
          gt_str_get(arguments->distlen_filename));
      had_err = -1;
    }
  }

  if (had_err != 0)
    gt_file_delete(arguments->outfp);

  return had_err;
}

static int gt_simreads_plot_disc_distri(unsigned long key,
                                        unsigned long long value,
                                        GtFile *outfile)
{
  gt_file_xprintf(outfile, "%lu %llu\n", key, value);
  return 0;
}

static int gt_simreads_write_dist_file(const char *title, GtDiscDistri *dist,
                                       GtStr *filename, GtError *err)
{
  int had_err = 0;
  GtFile *dfile;

  gt_assert(filename);
  gt_assert(dist);
  gt_error_check(err);
  dfile = gt_file_new(gt_str_get(filename), "w", err);
  if (!dfile)
  {
    had_err = -1;
  }
  else
  {
    gt_file_xprintf(dfile, "%s", title);
    gt_disc_distri_foreach(dist,
     (GtDiscDistriIterFunc) gt_simreads_plot_disc_distri, dfile);
    gt_file_delete(dfile);
  }
  return had_err;
}

typedef struct {
  unsigned long value;
  unsigned long length;
} GtSimreadsDistvalue;

static unsigned long gt_simreads_readlen_from_dist(
    const GtSimreadsDistvalue *dist, const unsigned long value,
    const unsigned long nofvalues)
{
  unsigned long l = 0, r = nofvalues - 1, m = r >> 1;
  while (value != dist[m].value)
  {
    if (value < dist[m].value)
    {
      gt_assert(l <= r);
      if (m == 0 || value > dist[m - 1].value)
        return dist[m].length;
      else
        r = m - 1;
    }
    else
    {
      gt_assert(m < nofvalues - 1);
      l = m + 1;
    }
    m = l + ((r - l) >> 1);
  }
  return dist[m].length;
}

static int gt_simreads_runner(GT_UNUSED int argc,
                              const char **argv, int parsed_args,
                              void *tool_arguments, GtError *err)
{
  GtSimreadsArguments *arguments = tool_arguments;
  GtEncseq *target = NULL;
  GtEncseqLoader *target_loader = NULL;
  GtEncseqReader *target_reader = NULL;
  GtReadmode readmode;
  GtStr *read = NULL, *description = NULL;
  unsigned long output_bases = 0, output_reads = 0,
                output_rcmode_reads = 0,
                required_output_bases = 0,
                target_total_length,
                startpos, readlen, i;
  unsigned long long progress = 0;
  bool fixed_readlen;
  int had_err = 0;
  GtUchar ch;
  GtLogger *logger = NULL;
  GtDiscDistri *lengths = NULL, *starts = NULL;
  GtSimreadsDistvalue *input_distlen = NULL;
  off_t j, nofvalues = 0;

  gt_error_check(err);
  gt_assert(arguments);

  logger = gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX, stderr);
  fixed_readlen = (arguments->maxlen == GT_UNDEF_ULONG &&
      arguments->minlen != GT_UNDEF_ULONG);
  if (arguments->ds)
    starts = gt_disc_distri_new();
  if (arguments->dl)
    lengths = gt_disc_distri_new();
  readlen = arguments->minlen;

  /* load encseq */
  target_loader = gt_encseq_loader_new();
  target = gt_encseq_loader_load(target_loader, argv[parsed_args], err);
  if (target == NULL)
    had_err = -1;
  if (!had_err) {
    target_reader = gt_encseq_create_reader_with_readmode(target,
                                                          GT_READMODE_FORWARD,
                                                          0);
  }

  if (!had_err)
  {
    target_total_length = gt_encseq_total_length(target);
    gt_logger_log(logger, "number of templates: %lu",
                  gt_encseq_num_of_sequences(target));
    gt_logger_log(logger, "total template length: %lu",
                  target_total_length);
    if (arguments->coverage != GT_UNDEF_ULONG)
    {
      gt_logger_log(logger, "required coverage: %lu",
                            arguments->coverage);
      required_output_bases = arguments->coverage * target_total_length;
      if (arguments->show_progressbar)
        gt_progressbar_start(&progress,
            (unsigned long long)required_output_bases);
    }
    else
    {
      gt_assert(arguments->num != GT_UNDEF_ULONG);
      gt_logger_log(logger, "required number of reads: %lu",
                            arguments->num);
      if (arguments->show_progressbar)
        gt_progressbar_start(&progress, (unsigned long long)arguments->num);
    }
    read = gt_str_new();
    description = gt_str_new();
    if (arguments->distlen_file != NULL)
    {
      gt_logger_log(logger, "read length distribution file: %s",
          gt_str_get(arguments->distlen_filename));
      nofvalues = gt_file_size(gt_str_get(arguments->distlen_filename)) /
                       (off_t) (sizeof (unsigned long) << 1);
      input_distlen = gt_malloc(sizeof(GtSimreadsDistvalue) * nofvalues);
      for (j = 0; j < nofvalues; j++)
      {
        (void)gt_xfread(&(input_distlen[j].length), sizeof (unsigned long),
            (size_t)1, arguments->distlen_file);
        (void)gt_xfread(&(input_distlen[j].value), sizeof (unsigned long),
            (size_t)1, arguments->distlen_file);
      }
      for (j = (off_t)1; j < nofvalues; j++)
      {
        input_distlen[j].value += input_distlen[j - 1].value;
      }
    }
    else if (fixed_readlen)
      gt_logger_log(logger, "required read length (fixed): %lu",
                    arguments->minlen);
    else
      gt_logger_log(logger, "required read length range: %lu-%lu",
                    arguments->minlen, arguments->maxlen);
  }

  while (!had_err)
  {
    gt_str_reset(read);
    if (!fixed_readlen)
    {
      if (arguments->distlen_file != NULL)
      {
        gt_assert(input_distlen != NULL);
        readlen = gt_simreads_readlen_from_dist(input_distlen,
            gt_rand_max(input_distlen[nofvalues - 1].value),
            (unsigned long)nofvalues);
        gt_assert(readlen <= input_distlen[nofvalues - 1].value);
      }
      else
      {
        readlen = gt_rand_max(arguments->maxlen -
            arguments->minlen) + arguments->minlen;
      }
      if (arguments->dl)
        gt_disc_distri_add(lengths, readlen);
    }
    gt_assert(target_total_length > readlen);
    startpos = gt_rand_max(target_total_length - readlen);
    readmode = (arguments->singlestrand || gt_rand_max(1UL)
      ? GT_READMODE_FORWARD : GT_READMODE_REVCOMPL);
    gt_encseq_reader_reinit_with_readmode(target_reader, target,
                                          readmode, startpos);
    for (i = 0; i < readlen; i++)
    {
      ch = gt_encseq_reader_next_encoded_char(target_reader);
      if (ch == (GtUchar)SEPARATOR)
      {
        break;
      }
      else
      {
        gt_str_append_char(read, gt_alphabet_decode(
                             gt_encseq_alphabet(target), ch));
      }
    }
    if (i < readlen)
    {
      /* a separator was found, discard read and restart */
      continue;
    }
    else
    {
      gt_str_reset(description);
      gt_str_append_cstr(description, "read_");
      gt_str_append_ulong(description, output_reads);
      gt_fasta_show_entry(gt_str_get(description),
                          gt_str_get(read),
                          gt_str_length(read),
                          60UL,
                          arguments->outfp);
      output_bases += gt_str_length(read);
      output_reads++;
      if (arguments->verbose && readmode == GT_READMODE_FORWARD)
        output_rcmode_reads++;
      if (arguments->ds)
      {
        if (readmode == GT_READMODE_FORWARD)
          gt_disc_distri_add(starts, startpos);
        else
          gt_disc_distri_add(starts, target_total_length - 1 - startpos);
      }
    }

    /* test break conditions and update progressbar */
    if (arguments->coverage != GT_UNDEF_ULONG)
    {
      if (arguments->show_progressbar)
        progress = (unsigned long long)output_bases;
      if (output_bases >= required_output_bases)
        break;
    }
    else
    {
      gt_assert(arguments->num != GT_UNDEF_ULONG);
      if (arguments->show_progressbar)
        progress = (unsigned long long)output_reads;
      if (output_reads == arguments->num)
        break;
    }
  }

  if (!had_err)
  {
    gt_logger_log(logger, "coverage: %.3f",
                  (float) output_bases / target_total_length);
    gt_logger_log(logger, "total reads length: %lu", output_bases);
    if (!fixed_readlen)
      gt_logger_log(logger, "average reads length: %.1f",
                    (float) output_bases / output_reads);
    gt_logger_log(logger, "number of reads: %lu", output_reads);
    gt_logger_log(logger, "- forward: %lu", output_reads-output_rcmode_reads);
    gt_logger_log(logger, "- revcompl: %lu", output_rcmode_reads);
  }

  if (!had_err && arguments->dl)
    had_err = gt_simreads_write_dist_file(
                   "# distribution of read lengths:\n", lengths,
                   arguments->dlfilename, err);
  if (!had_err && arguments->ds)
    had_err = gt_simreads_write_dist_file(
                   "# distribution of start positions:\n", starts,
                   arguments->dsfilename, err);

  if (arguments->show_progressbar)
    gt_progressbar_stop();
  gt_free(input_distlen);
  gt_str_delete(read);
  gt_str_delete(description);
  gt_encseq_reader_delete(target_reader);
  gt_encseq_loader_delete(target_loader);
  gt_encseq_delete(target);
  gt_logger_delete(logger);
  if (arguments->dl)
    gt_disc_distri_delete(lengths);
  if (arguments->ds)
    gt_disc_distri_delete(starts);
  return had_err;
}

GtTool* gt_simreads(void)
{
  return gt_tool_new(gt_simreads_arguments_new,
                     gt_simreads_arguments_delete,
                     gt_simreads_option_parser_new,
                     gt_simreads_arguments_check,
                     gt_simreads_runner);
}
