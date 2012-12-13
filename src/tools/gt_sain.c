/*
  Copyright (c) 2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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
#include "core/unused_api.h"
#include "core/fa.h"
#include "core/timer_api.h"
#include "core/showtime.h"
#include "core/logger.h"
#include "tools/gt_sain.h"
#include "match/sfx-sain.h"

typedef struct
{
  bool icheck, fcheck, verbose;
  GtStr *encseqfile, *plainseqfile, *dir;
  GtReadmode readmode;
} GtSainArguments;

static void* gt_sain_arguments_new(void)
{
  GtSainArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->encseqfile = gt_str_new();
  arguments->plainseqfile = gt_str_new();
  arguments->dir = gt_str_new_cstr("fwd");
  arguments->readmode = GT_READMODE_FORWARD;
  return arguments;
}

static void gt_sain_arguments_delete(void *tool_arguments)
{
  GtSainArguments *arguments = tool_arguments;

  if (arguments != NULL)
  {
    gt_str_delete(arguments->encseqfile);
    gt_str_delete(arguments->plainseqfile);
    gt_str_delete(arguments->dir);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_sain_option_parser_new(void *tool_arguments)
{
  GtSainArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *optionfcheck, *optionesq, *optionfile;

  gt_assert(arguments != NULL);

  /* init */
  op = gt_option_parser_new("[option ...] [file]", /* XXX */
                            "Compute suffix array using induced "
                            "suffix sorting.");
  gt_option_parser_set_mail_address(op,"<kurtz@zbh.uni-hamburg.de>");

  /* -esq */
  optionesq = gt_option_new_string("esq", "specify encseq file",
                             arguments->encseqfile, NULL);
  gt_option_parser_add_option(op, optionesq);

  /* -dir */
  gt_encseq_options_add_readmode_option(op, arguments->dir);

  /* -file */
  optionfile = gt_option_new_string("file", "specify filename",
                                   arguments->plainseqfile, NULL);
  gt_option_parser_add_option(op, optionfile);

  /* -icheck */
  option = gt_option_new_bool("icheck",
                              "intermediate check of all sorted arrays",
                              &arguments->icheck, false);
  gt_option_parser_add_option(op, option);

  /* -fcheck */
  optionfcheck = gt_option_new_bool("fcheck", "final check of suffix array",
                                    &arguments->fcheck, false);
  gt_option_parser_add_option(op, optionfcheck);
  gt_option_imply(optionfcheck, optionesq);
  gt_option_exclude(optionesq,optionfile);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_sain_option_parser_check(int rest_argc,
                                       void *tool_arguments, GtError *err)
{
  int retval;
  GtSainArguments *arguments = tool_arguments;

  if (rest_argc > 0)
  {
    gt_error_set(err,"%d superfluous argument%s",
                 rest_argc,rest_argc > 1 ? "s" : "");
    return -1;
  }
  retval = gt_readmode_parse(gt_str_get(arguments->dir), err);
  if (retval < 0)
  {
    return -1;
  } else
  {
    arguments->readmode = (GtReadmode) retval;
    return 0;
  }
}

static int gt_sain_checkmaxsequencelength(unsigned long len,bool forencseq,
                                          GtError *err)
{
  unsigned long maxsequencelength;

  if (forencseq)
  {
    maxsequencelength = (unsigned long) (~GT_FIRSTBIT) - 1 - GT_COMPAREOFFSET;
  } else
  {
    maxsequencelength = (unsigned long) (~GT_FIRSTBIT) - 1;
  }
  if (len > maxsequencelength)
  {
    gt_error_set(err,"sequence of size %lu is too long: sain algorithm "
                     "can only compute sequence of length up to %lu",
                     len,maxsequencelength);
    return -1;
  }
  return 0;
}

typedef struct
{
  GtTimer *timer;
  GtLogger *logger;
} GtSainTimerandLogger;

static GtSainTimerandLogger *gt_sain_timer_logger_new(bool verbose)
{
  GtSainTimerandLogger *tl = gt_malloc(sizeof (*tl));

  tl->timer = NULL;
  tl->logger = gt_logger_new(verbose,GT_LOGGER_DEFLT_PREFIX,stdout);
  if (gt_showtime_enabled())
  {
    if (verbose)
    {
      tl->timer = gt_timer_new_with_progress_description(
                                       "allocate suftab and undef entries");
    } else
    {
      tl->timer = gt_timer_new();
      gt_timer_omit_last_stage(tl->timer);
    }
    gt_timer_start(tl->timer);
  }
  return tl;
}

static void gt_sain_timer_logger_delete(GtSainTimerandLogger *tl)
{
  if (tl->timer != NULL)
  {
    gt_timer_show_progress_final(tl->timer, stdout);
    gt_timer_stop(tl->timer);
    gt_timer_delete(tl->timer);
  }
  gt_logger_delete(tl->logger);
  gt_free(tl);
}

static int gt_sain_runner(int argc, GT_UNUSED const char **argv,
                          int parsed_args, void *tool_arguments, GtError *err)
{
  GtSainArguments *arguments = tool_arguments;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments != NULL);

  gt_assert (argc >= parsed_args);
  if (parsed_args < argc)
  {
    gt_error_set(err,"superfluous arguments");
    had_err = -1;
  } else
  {
    if (gt_str_length(arguments->encseqfile) > 0)
    {
      GtEncseqLoader *el = gt_encseq_loader_new();
      GtEncseq *encseq
        = gt_encseq_loader_load(el,gt_str_get(arguments->encseqfile),err);
      if (encseq == NULL)
      {
        had_err = -1;
      } else
      {
        if (gt_sain_checkmaxsequencelength(gt_encseq_total_length(encseq),true,
                                           err)
                       != 0)
        {
          had_err = -1;
        } else
        {
          if (!gt_alphabet_is_dna(gt_encseq_alphabet(encseq)) &&
              (arguments->readmode == GT_READMODE_COMPL ||
               arguments->readmode == GT_READMODE_REVCOMPL))
          {
            gt_error_set(err,"option -dir cpl and -dir rcl is only "
                             "possible for DNA sequences");
            had_err = -1;
          }
        }
        if (!had_err)
        {
          GtSainTimerandLogger *tl
            = gt_sain_timer_logger_new(arguments->verbose);
          gt_sain_encseq_sortsuffixes(encseq,
                                      arguments->readmode,
                                      arguments->icheck,
                                      arguments->fcheck,
                                      tl->logger,
                                      tl->timer);
          gt_sain_timer_logger_delete(tl);
        }
      }
      gt_encseq_delete(encseq);
      gt_encseq_loader_delete(el);
    }
    if (gt_str_length(arguments->plainseqfile) > 0)
    {
      GtUchar *plainseq;
      size_t len;

      plainseq = gt_fa_mmap_read(gt_str_get(arguments->plainseqfile),&len,err);
      if (plainseq == NULL)
      {
        had_err = -1;
      } else
      {
        if (gt_sain_checkmaxsequencelength((unsigned long) len,false,err) != 0)
        {
          had_err = -1;
        } else
        {
          if (arguments->readmode != GT_READMODE_FORWARD)
          {
            gt_error_set(err,"option -dir and -file exclude each other");
            had_err = -1;
          }
        }
        if (!had_err)
        {
          GtSainTimerandLogger *tl
            = gt_sain_timer_logger_new(arguments->verbose);
          gt_sain_plain_sortsuffixes(plainseq,
                                     (unsigned long) len,
                                     arguments->icheck,
                                     tl->logger,
                                     tl->timer);
          gt_sain_timer_logger_delete(tl);
        }
      }
      gt_fa_xmunmap(plainseq);
    }
  }
  return had_err;
}

GtTool* gt_sain(void)
{
  return gt_tool_new(gt_sain_arguments_new,
                     gt_sain_arguments_delete,
                     gt_sain_option_parser_new,
                     gt_sain_option_parser_check,
                     gt_sain_runner);
}
