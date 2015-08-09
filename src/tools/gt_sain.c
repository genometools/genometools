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
#include "core/xansi_api.h"
#include "core/basename_api.h"
#include "tools/gt_sain.h"
#include "match/sfx-linlcp.h"
#include "match/sfx-sain.h"
#include "match/bare-encseq.h"

typedef struct
{
  bool icheck, fcheck, outlcptab, lcpkasai,
       verbose, dommap, dnaalphabet, proteinalphabet, suftabout, tistabout;
  GtStr *encseqfile, *plainseqfile, *fastafile, *dir, *smap;
  GtReadmode readmode;
} GtSainArguments;

static void* gt_sain_arguments_new(void)
{
  GtSainArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->suftabout = false;
  arguments->tistabout = false;
  arguments->encseqfile = gt_str_new();
  arguments->plainseqfile = gt_str_new();
  arguments->fastafile = gt_str_new();
  arguments->smap = gt_str_new();
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
    gt_str_delete(arguments->fastafile);
    gt_str_delete(arguments->smap);
    gt_str_delete(arguments->dir);
    gt_free(arguments);
  }
}

static GtOptionParser *gt_sain_option_parser_new(void *tool_arguments)
{
  GtSainArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *optionesq, *optionfile, *optionmmap,
           *optionfasta, *optiondnaalphabet, *optionproteinalphabet,
           *optionlcp, *optionlcpkasai, *optionsmap, *optiontistabout;

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

  /* -lcp */
  optionlcp = gt_option_new_bool("lcp", "output lcp table",
                                 &arguments->outlcptab, false);
  gt_option_parser_add_option(op, optionlcp);

  /* -kasai */
  optionlcpkasai = gt_option_new_bool("kasai", "use kasai algorithm to compute "
                                      "lcp table",
                                      &arguments->lcpkasai, false);
  gt_option_parser_add_option(op, optionlcpkasai);

  /* -file */
  optionfile = gt_option_new_string("file", "specify filename",
                                   arguments->plainseqfile, NULL);
  gt_option_parser_add_option(op, optionfile);

  /* -fasta */
  optionfasta = gt_option_new_string("fasta","fasta input",
                                    arguments->fastafile, NULL);
  gt_option_parser_add_option(op, optionfasta);

  /* -dna */
  optiondnaalphabet = gt_option_new_bool("dna","use DNA alphabet",
                                         &arguments->dnaalphabet, false);
  gt_option_parser_add_option(op, optiondnaalphabet);

  /* -protein */
  optionproteinalphabet = gt_option_new_bool("protein","use protein alphabet",
                                             &arguments->proteinalphabet,false);
  gt_option_parser_add_option(op, optionproteinalphabet);

  /* -smap */
  optionsmap = gt_option_new_string("smap","specify symbol map file",
                                    arguments->smap, NULL);
  gt_option_parser_add_option(op, optionsmap);

  /* -icheck */
  option = gt_option_new_bool("icheck",
                              "intermediate check of all sorted arrays",
                              &arguments->icheck, false);
  gt_option_parser_add_option(op, option);

  /* -mmap */
  optionmmap = gt_option_new_bool("mmap",
                                  "use mmap to map the input file "
                                  "(requires option -file)",
                                  &arguments->dommap, false);
  gt_option_parser_add_option(op, optionmmap);

  /* -fcheck */
  option = gt_option_new_bool("fcheck", "final check of suffix array",
                              &arguments->fcheck, false);
  gt_option_parser_add_option(op, option);

  /* -suf */
  option = gt_option_new_bool("suf","output suffix array (table suftab)",
                              &arguments->suftabout, false);
  gt_option_parser_add_option(op, option);

  /* -tis */
  optiontistabout
    = gt_option_new_bool("tis","output transformed input sequence",
                         &arguments->tistabout, false);
  gt_option_parser_add_option(op, optiontistabout);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  gt_option_exclude(optionesq,optionfile);
  gt_option_exclude(optionesq, optionfasta);
  gt_option_exclude(optionfile, optionfasta);
  gt_option_exclude(optionfasta, optionmmap);
  gt_option_exclude(optiondnaalphabet,optionproteinalphabet);
  gt_option_exclude(optiondnaalphabet,optionsmap);
  gt_option_exclude(optionproteinalphabet,optionsmap);
  gt_option_imply(optiondnaalphabet,optionfasta);
  gt_option_imply(optionlcpkasai,optionlcp);
  gt_option_imply(optionproteinalphabet,optionfasta);
  gt_option_imply(optionsmap,optionfasta);
  gt_option_imply_either_2(optiontistabout,optionfasta,optionfile);
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

static GtAlphabet *gt_sain_getalphabet(const GtSainArguments *arguments,
                                       GtError *err)
{
  GtAlphabet *alphabet;
  int had_err = 0;

  if (arguments->dnaalphabet)
  {
    alphabet = gt_alphabet_new_dna();
  } else
  {
    if (arguments->proteinalphabet)
    {
      alphabet = gt_alphabet_new_protein();
    } else
    {
      if (gt_str_length(arguments->smap) > 0)
      {
        alphabet = gt_alphabet_new_from_file_no_suffix(
                                           gt_str_get(arguments->smap),
                                           err);
      } else
      {
        GtStrArray *filenametab = gt_str_array_new();

        gt_assert(gt_str_length(arguments->fastafile) > 0);
        gt_str_array_add_cstr(filenametab,gt_str_get(arguments->fastafile));
        alphabet = gt_alphabet_new_from_sequence(filenametab, err);
        gt_str_array_delete(filenametab);
      }
      if (alphabet == NULL)
      {
        had_err = -1;
      } else
      {
        if (!gt_alphabet_is_dna(alphabet) &&
            (arguments->readmode == GT_READMODE_COMPL ||
             arguments->readmode == GT_READMODE_REVCOMPL))
        {
          gt_error_set(err,"option -dir cpl and -dir rcl is only "
                           "possible for DNA sequences");
          had_err = -1;
        }
      }
    }
  }
  if (had_err)
  {
    gt_alphabet_delete(alphabet);
    return NULL;
  }
  gt_assert(alphabet != NULL);
  return alphabet;
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
    GtUsainindextype *suftab = NULL;

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
          suftab = gt_sain_encseq_sortsuffixes(encseq,
                                               arguments->readmode,
                                               arguments->icheck,
                                               arguments->fcheck,
                                               tl->logger,
                                               tl->timer);
          gt_sain_timer_logger_delete(tl);
          gt_free(suftab);
        }
      }
      gt_encseq_delete(encseq);
      gt_encseq_loader_delete(el);
    }
    if (!had_err && (gt_str_length(arguments->plainseqfile) > 0 ||
                     gt_str_length(arguments->fastafile) > 0))
    {
      GtUchar *filecontents;
      size_t len;
      const char *inputfile = gt_str_length(arguments->plainseqfile) > 0
                                ? gt_str_get(arguments->plainseqfile)
                                : gt_str_get(arguments->fastafile);

      if (arguments->dommap)
      {
        filecontents = gt_fa_mmap_read (inputfile,&len,err);
      } else
      {
        filecontents = gt_fa_heap_read (inputfile,&len,err);
      }
      if (filecontents == NULL)
      {
        had_err = -1;
      } else
      {
        GtBareEncseq *bare_encseq = NULL;

        if (!had_err)
        {
          if (gt_str_length(arguments->fastafile) > 0)
          {
            GtAlphabet *alphabet = gt_sain_getalphabet(arguments,err);

            if (alphabet == NULL)
            {
              had_err = -1;
            } else
            {
              bare_encseq = gt_bare_encseq_parse_new(filecontents,len,alphabet,
                                                     err);
              if (bare_encseq == NULL)
              {
                had_err = -1;
              }
            }
            gt_alphabet_delete(alphabet);
          }
        }
        if (!had_err &&
            gt_sain_checkmaxsequencelength((GtUword) len,false,err) != 0)
        {
          had_err = -1;
        }
        if (!had_err)
        {
          GtSainTimerandLogger *tl
            = gt_sain_timer_logger_new(arguments->verbose);
          if (gt_str_length(arguments->plainseqfile) > 0)
          {
            suftab = gt_sain_plain_sortsuffixes(filecontents,
                                                (GtUword) len,
                                                arguments->icheck,
                                                arguments->fcheck,
                                                tl->logger,
                                                tl->timer);
          } else
          {
            gt_assert(gt_str_length(arguments->fastafile) > 0 &&
                      bare_encseq != NULL);
            if (arguments->readmode != GT_READMODE_FORWARD)
            {
              bare_encseq_convert(bare_encseq,
                                  GT_ISDIRREVERSE(arguments->readmode)
                                    ? false : true,
                                  GT_ISDIRCOMPLEMENT(arguments->readmode)
                                    ? false : true);
            }
            suftab = gt_sain_bare_encseq_sortsuffixes(bare_encseq,
                                                      arguments->readmode,
                                                      arguments->icheck,
                                                      arguments->fcheck,
                                                      tl->logger,
                                                      tl->timer);
          }
          if (arguments->suftabout)
          {
            char *inputfile_basename = gt_basename(inputfile);
            FILE *fpout = gt_fa_fopen_with_suffix(inputfile_basename,
                                                  ".suf","wb",err);
            if (fpout == NULL)
            {
              had_err = -1;
            } else
            {
              GtUword totallength;
              if (bare_encseq != NULL)
              {
                totallength = gt_bare_encseq_total_length(bare_encseq);
              } else
              {
                totallength = len;
              }
              gt_xfwrite(suftab,sizeof *suftab,totallength+1,fpout);
              gt_fa_fclose(fpout);
            }
            gt_free(inputfile_basename);
          }
          if (arguments->tistabout)
          {
            char *inputfile_basename = gt_basename(inputfile);
            FILE *fpout = gt_fa_fopen_with_suffix(inputfile_basename,
                                                  ".tis","wb",err);
            if (fpout == NULL)
            {
              had_err = -1;
            } else
            {
              GtUword totallength;
              const GtUchar *sequence;

              if (bare_encseq != NULL)
              {
                totallength = gt_bare_encseq_total_length(bare_encseq);
                sequence = gt_bare_encseq_sequence(bare_encseq);
              } else
              {
                totallength = len;
                sequence = filecontents;
              }
              gt_xfwrite(sequence,sizeof *sequence,totallength,fpout);
              gt_fa_fclose(fpout);
            }
            gt_free(inputfile_basename);
          }
          if (!had_err && arguments->outlcptab)
          {
            unsigned int *lcptab;
            GtUword maxlcp = 0, partwidth, totallength;
            bool withspecial;

            if (gt_str_length(arguments->plainseqfile) > 0)
            {
              withspecial = false;
              totallength = (GtUword) len;
              partwidth = totallength;
            } else
            {
              withspecial = true;
              totallength = gt_bare_encseq_total_length(bare_encseq);
              partwidth = totallength -
                          gt_bare_encseq_specialcharacters(bare_encseq);
            }
            if (tl->timer != NULL && gt_logger_enabled(tl->logger))
            {
              gt_timer_show_progress(tl->timer,"compute lcp-value",stdout);
            }
            if (arguments->lcpkasai)
            {
              lcptab = gt_plain_lcp13_kasai(&maxlcp,
                                            filecontents,
                                            withspecial,
                                            partwidth,
                                            totallength,
                                            suftab);
            } else
            {
              lcptab = gt_plain_lcp_phialgorithm(false,
                                                 &maxlcp,
                                                 filecontents,
                                                 withspecial,
                                                 partwidth,
                                                 totallength,
                                                 suftab);
            }
            gt_logger_log(tl->logger,"maxlcp="GT_WU,maxlcp);
            gt_free(lcptab);
          }
          gt_free(suftab);
          gt_sain_timer_logger_delete(tl);
        }
        gt_bare_encseq_delete(bare_encseq);
      }
      if (arguments->dommap)
      {
        gt_fa_xmunmap(filecontents);
      } else
      {
        gt_free(filecontents);
      }
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
