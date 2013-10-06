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

#ifndef S_SPLINT_S
#include <ctype.h>
#endif
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
  GtUchar *sequence;
  GtUword totallength, specialcharacters, numofchars, *charcount;
} GtBareEncseq;

static void gt_bare_encseq_delete(GtBareEncseq *bare_encseq)
{
  if (bare_encseq != NULL)
  {
    gt_free(bare_encseq->charcount);
    gt_free(bare_encseq);
  }
}

static GtBareEncseq *gt_bare_encseq_new(GtUchar *filecontents,
                                        size_t numofbytes,
                                        GtError *err)
{
  GtUchar *writeptr = filecontents, *readptr = filecontents,
          smap[UCHAR_MAX+1];
  size_t idx;
  const GtUchar undefined = (GtUchar) UCHAR_MAX,
        *endptr = filecontents + numofbytes;
  bool firstline = true, haserr = false;
  const char *wildcard_list = "nsywrkvbdhmNSYWRKVBDHM";
  GtBareEncseq *bare_encseq = gt_malloc(sizeof *bare_encseq);

  for (idx = 0; idx <= (size_t) UCHAR_MAX; idx++)
  {
    smap[idx] = undefined;
  }
  smap['a'] = 0;
  smap['A'] = 0;
  smap['c'] = (GtUchar) 1;
  smap['C'] = (GtUchar) 1;
  smap['g'] = (GtUchar) 2;
  smap['G'] = (GtUchar) 2;
  smap['t'] = (GtUchar) 3;
  smap['T'] = (GtUchar) 3;
  smap['u'] = (GtUchar) 3;
  smap['U'] = (GtUchar) 3;
  bare_encseq->specialcharacters = 0;
  bare_encseq->numofchars = 4UL;
  bare_encseq->charcount = gt_calloc((size_t) bare_encseq->numofchars,
                                     sizeof *bare_encseq->charcount);
  for (idx = 0; idx < strlen(wildcard_list); idx++)
  {
    smap[(int) wildcard_list[idx]] = WILDCARD;
  }
  readptr = filecontents;
  while (!haserr && readptr < endptr)
  {
    if (*readptr == '>')
    {
      if (!firstline)
      {
        *writeptr++ = SEPARATOR;
        bare_encseq->specialcharacters++;
      } else
      {
        firstline = false;
      }
      while (readptr < endptr && *readptr != '\n')
      {
        readptr++;
      }
      readptr++;
    } else
    {
      while (readptr < endptr && *readptr != '\n')
      {
        if (!isspace(*readptr))
        {
          GtUchar cc = smap[*readptr];
          if (cc == undefined)
          {
            gt_error_set(err,"illegal input characters %c\n",*readptr);
            haserr = true;
            break;
          }
          if (ISSPECIAL(cc))
          {
            bare_encseq->specialcharacters++;
          } else
          {
            gt_assert((GtUword) cc < bare_encseq->numofchars);
            bare_encseq->charcount[(int) cc]++;
          }
          *writeptr++ = cc;
        }
        readptr++;
      }
      readptr++;
    }
  }
  bare_encseq->sequence = filecontents;
  bare_encseq->totallength = (GtUword) (writeptr - filecontents);
  if (haserr)
  {
    gt_bare_encseq_delete(bare_encseq);
    return NULL;
  }
  return bare_encseq;
}

typedef struct
{
  bool icheck, fcheck, verbose, dommap;
  GtStr *encseqfile, *plainseqfile, *fastadnafile, *dir;
  GtReadmode readmode;
} GtSainArguments;

static void* gt_sain_arguments_new(void)
{
  GtSainArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->encseqfile = gt_str_new();
  arguments->plainseqfile = gt_str_new();
  arguments->fastadnafile = gt_str_new();
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
    gt_str_delete(arguments->fastadnafile);
    gt_str_delete(arguments->dir);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_sain_option_parser_new(void *tool_arguments)
{
  GtSainArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *optionfcheck, *optionesq, *optionfile, *optionmmap,
           *optionfastadna;

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

  /* -fastadna */
  optionfastadna = gt_option_new_string("fastadna",
                                        "fasta input with DNA sequence",
                                        arguments->fastadnafile, NULL);
  gt_option_parser_add_option(op, optionfastadna);

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
  optionfcheck = gt_option_new_bool("fcheck", "final check of suffix array",
                                    &arguments->fcheck, false);
  gt_option_parser_add_option(op, optionfcheck);
  gt_option_imply(optionfcheck, optionesq);
  gt_option_exclude(optionesq,optionfile);
  gt_option_exclude(optionesq, optionfastadna);
  gt_option_exclude(optionfile, optionfastadna);
  gt_option_exclude(optionfastadna, optionmmap);
  gt_option_imply(optionmmap, optionfile);

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

static int gt_sain_checkmaxsequencelength(GtUword len,bool forencseq,
                                          GtError *err)
{
  GtUword maxsequencelength;

  if (forencseq)
  {
    maxsequencelength = (GtUword) (~GT_FIRSTBIT) - 1 - GT_COMPAREOFFSET;
  } else
  {
    maxsequencelength = (GtUword) (~GT_FIRSTBIT) - 1;
  }
  if (len > maxsequencelength)
  {
    gt_error_set(err,"sequence of size "GT_WU" is too long: sain algorithm "
                     "can only compute sequence of length up to "GT_WU"",
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
    if (gt_str_length(arguments->plainseqfile) > 0 ||
        gt_str_length(arguments->fastadnafile) > 0)
    {
      GtUchar *filecontents;
      size_t len;
      const char *inputfile = gt_str_length(arguments->plainseqfile) > 0
                                ? gt_str_get(arguments->plainseqfile)
                                : gt_str_get(arguments->fastadnafile);

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

        if (gt_str_length(arguments->fastadnafile) > 0)
        {
          bare_encseq = gt_bare_encseq_new(filecontents,len,err);
          if (bare_encseq == NULL)
          {
            had_err = -1;
          }
        }
        if (!had_err)
        {
          if (gt_sain_checkmaxsequencelength((GtUword) len,false,err) != 0)
          {
            had_err = -1;
          } else
          {
            if (arguments->readmode != GT_READMODE_FORWARD)
            {
              gt_error_set(err,"option -dir and -file/-fastadna exclude "
                               " each other");
              had_err = -1;
            }
          }
        }
        if (!had_err)
        {
          GtSainTimerandLogger *tl
            = gt_sain_timer_logger_new(arguments->verbose);
          if (gt_str_length(arguments->plainseqfile) > 0)
          {
            gt_sain_plain_sortsuffixes(filecontents,
                                       (GtUword) len,
                                       arguments->icheck,
                                       tl->logger,
                                       tl->timer);
          }
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
