/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2010-2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2007-2011 Center for Bioinformatics, University of Hamburg

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
#include <limits.h>
#ifndef S_SPLINT_S
#include <sys/types.h>
#endif
#include "core/fileutils_api.h"
#include "core/xposix.h"
#include "core/fa.h"
#include "core/ma.h"
#include "core/versionfunc.h"
#include "core/seq_iterator_sequence_buffer_api.h"
#include "core/seq_iterator_fastq_api.h"
#include "core/unused_api.h"
#include "core/versionfunc.h"
#include "core/xansi_api.h"
#include "core/progressbar.h"
#include "core/format64.h"
#include "core/unused_api.h"
#include "core/types_api.h"
#include "core/disc_distri_api.h"
#include "core/option_api.h"
#include "core/error_api.h"
#include "core/logger.h"
#include "extended/assembly_stats_calculator.h"
#include "match/stamp.h"
#include "tools/gt_seqstat.h"

typedef struct
{
  bool verbose,
       dodistlen,
       binarydistlen,
       doastretch,
       docstats,
       showestimsize;
  unsigned int bucketsize;
  unsigned long genome_length;
} SeqstatArguments;

#define GT_SEQSTAT_BINARY_DISTLEN_SUFFIX ".distlen"

static GtOPrval parse_options(SeqstatArguments *arguments,
                              int *parsed_args,int argc,
                              const char **argv, GtError *err)
{
  GtOptionParser *op;
  GtOption *optionverbose, *optiondistlen, *optionbucketsize,
           *optioncontigs, *optionastretch, *optionestimsize,
           *optionbinarydistlen, *optiongenome;
  GtOPrval oprval;

  gt_error_check(err);

  op = gt_option_parser_new("[options] file [...]",
                            "Calculate statistics for fasta file(s).");

  optionverbose = gt_option_new_bool("v","be verbose",
                                  &arguments->verbose,false);
  gt_option_parser_add_option(op, optionverbose);

  optiondistlen = gt_option_new_bool("distlen",
                                  "show distribution of sequence length",
                                  &arguments->dodistlen,false);
  gt_option_parser_add_option(op, optiondistlen);

  optionbucketsize = gt_option_new_uint_min("b",
                                "bucket size for distlen option",
                                &arguments->bucketsize, 100U, 1U);
  gt_option_imply(optionbucketsize, optiondistlen);
  gt_option_parser_add_option(op, optionbucketsize);

  optionbinarydistlen = gt_option_new_bool("binary",
                                  "use a binary format for distlen output\n"
                                  "output filename: <first_input_filename>"
                                  GT_SEQSTAT_BINARY_DISTLEN_SUFFIX "\n"
                                  "bucketsize: 1",
                                  &arguments->binarydistlen,false);
  gt_option_imply(optionbinarydistlen, optiondistlen);
  gt_option_exclude(optionbinarydistlen, optionbucketsize);
  gt_option_parser_add_option(op, optionbinarydistlen);

  optioncontigs = gt_option_new_bool("contigs",
                                   "summary of contigs set statistics",
                                   &arguments->docstats,false);
  gt_option_parser_add_option(op, optioncontigs);

  optionastretch = gt_option_new_bool("astretch",
                                   "show distribution of A-substrings",
                                   &arguments->doastretch,false);
  gt_option_exclude(optiondistlen, optionastretch);
  gt_option_parser_add_option(op, optionastretch);
  gt_option_is_extended_option(optionastretch);

  optionestimsize = gt_option_new_bool("estimsize",
                                   "show estimated size",
                                   &arguments->showestimsize,false);
  gt_option_parser_add_option(op, optionestimsize);
  gt_option_is_development_option(optionestimsize);

  optiongenome = gt_option_new_ulong_min("genome",
                                "set genome length for NG50/NG80 calculation",
                                &arguments->genome_length, 0, 1UL);
  gt_option_imply(optiongenome, optioncontigs);
  gt_option_parser_add_option(op, optiongenome);

  gt_option_parser_set_min_args(op, 1U);
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gt_versionfunc,
                                  err);
  gt_option_parser_delete(op);
  return oprval;
}

static void showdistseqlen(unsigned long key, unsigned long long value,
                            void *data)
{
  unsigned long distvalue;
  unsigned int *bucketsize = data;

  gt_assert(value <= (unsigned long long) ULONG_MAX);
  distvalue = (unsigned long) value;
  printf("%lu--%lu %lu\n",
         (*bucketsize) * key,
         (*bucketsize) * (key+1) - 1,
         distvalue);
}

static void showdistseqlenbinary(unsigned long key, unsigned long long value,
                            void *data)
{
  unsigned long value_ul;
  /*@ignore@*/
  FILE *file = data;
  /*@end@*/

  gt_assert(value <= (unsigned long long) ULONG_MAX);
  value_ul = (unsigned long) value;
  gt_xfwrite(&key, sizeof (key), (size_t)1, file);
  gt_xfwrite(&value_ul, sizeof (value_ul), (size_t)1, file);
}

typedef struct
{
  unsigned long long sumA, *mmercount;
  unsigned long maxvalue, minkey;
} Astretchinfo;

static void showastretches(unsigned long key, unsigned long long value,
                           void *data)
{
  Astretchinfo *astretchinfo = (Astretchinfo *) data;

  astretchinfo->sumA += value * (unsigned long long) key;
  if (key > astretchinfo->maxvalue)
  {
    astretchinfo->maxvalue = key;
  }
  /*@ignore@*/
  printf("%lu %llu\n", key,value);
  /*@end@*/
}

static void showmeroccurrence(unsigned long key, unsigned long long value,
                              void *data)
{
  unsigned long len;
  Astretchinfo *astretchinfo = (Astretchinfo *) data;

  for (len=astretchinfo->minkey; len<= key; len++)
  {
    gt_assert(len <= astretchinfo->maxvalue);
    astretchinfo->mmercount[len] += value * (unsigned long long) (key-len+1);
  }
}

static unsigned long long accumulateastretch(GtDiscDistri *distastretch,
                                             const GtUchar *sequence,
                                             unsigned long len)
{
  unsigned long i, lenofastretch = 0;
  unsigned long long countA = 0;

  for (i=0; i<len; i++)
  {
    if (sequence[i] == 'A' || sequence[i] == 'a')
    {
      countA++;
      lenofastretch++;
    } else
    {
      if (lenofastretch > 0)
      {
        gt_disc_distri_add(distastretch,lenofastretch);
        lenofastretch = 0;
      }
    }
  }
  if (lenofastretch > 0)
  {
    gt_disc_distri_add(distastretch,lenofastretch);
  }
  return countA;
}

static void processastretches(const GtDiscDistri *distastretch,
                              GT_UNUSED unsigned long long countA)
{
  Astretchinfo astretchinfo;
  unsigned long len;

  astretchinfo.sumA = 0;
  astretchinfo.maxvalue = 0;
  astretchinfo.minkey = 10UL;
  gt_disc_distri_foreach(distastretch,showastretches,&astretchinfo);
  astretchinfo.mmercount = gt_malloc(sizeof (*astretchinfo.mmercount) *
                                    (astretchinfo.maxvalue+1));
  memset(astretchinfo.mmercount,0,sizeof (*astretchinfo.mmercount) *
                                  (astretchinfo.maxvalue+1));
  gt_disc_distri_foreach(distastretch,showmeroccurrence,&astretchinfo);
  for (len=astretchinfo.minkey; len<=astretchinfo.maxvalue; len++)
  {
    /*@ignore@*/
    printf("a^{%lu} occurs %llu times\n", len,astretchinfo.mmercount[len]);
    /*@end@*/
  }
  gt_assert(astretchinfo.sumA == countA);
  gt_free(astretchinfo.mmercount);
}

int gt_seqstat(int argc, const char **argv, GtError *err)
{
  GtStrArray *files;
  GtSeqIterator *seqit;
  const GtUchar *sequence;
  char *desc;
  unsigned long len;
  int i, parsed_args, had_err = 0;
  off_t totalsize;
  GtDiscDistri *distseqlen = NULL;
  GtAssemblyStatsCalculator *asc = NULL;
  GtLogger *asc_logger = NULL;
  GtDiscDistri *distastretch = NULL;
  uint64_t numofseq = 0, sumlength = 0;
  unsigned long minlength = 0, maxlength = 0;
  unsigned long long countA = 0;
  bool minlengthdefined = false;
  SeqstatArguments arguments;

  gt_error_check(err);

  /* option parsing */
  switch (parse_options(&arguments,&parsed_args, argc, argv, err)) {
    case GT_OPTION_PARSER_OK: break;
    case GT_OPTION_PARSER_ERROR:
        return -1;
    case GT_OPTION_PARSER_REQUESTS_EXIT:
        return 0;
  }

  files = gt_str_array_new();
  for (i = parsed_args; i < argc; i++)
  {
    gt_str_array_add_cstr(files, argv[i]);
  }
  totalsize = gt_files_estimate_total_size(files);
  if (arguments.showestimsize)
  {
    printf("# estimated total size is " Formatuint64_t "\n",
              PRINTuint64_tcast(totalsize));
  }
  if (!had_err) {
    /* read input using seqiterator */
    seqit = gt_seq_iterator_sequence_buffer_new(files, err);
    if (!seqit)
      had_err = -1;
    if (!had_err)
    {
      if (arguments.dodistlen)
      {
        distseqlen = gt_disc_distri_new();
        if (arguments.binarydistlen)
          arguments.bucketsize = 1U;
      }
      if (arguments.docstats)
      {
        asc = gt_assembly_stats_calculator_new();
        gt_assembly_stats_calculator_set_genome_length(asc,
            arguments.genome_length);
      }
      if (arguments.doastretch)
      {
        distastretch = gt_disc_distri_new();
      }
      if (arguments.verbose)
      {
        gt_progressbar_start(gt_seq_iterator_getcurrentcounter(seqit,
                                                            (unsigned long long)
                                                            totalsize),
                             (unsigned long long) totalsize);
      }
      while (true)
      {
        desc = NULL;
        had_err = gt_seq_iterator_next(seqit, &sequence, &len, &desc, err);
        if (had_err != 1) break; /* 0: finished; 1: error */
        if (arguments.dodistlen || arguments.docstats)
        {
          if (!minlengthdefined || minlength > len)
          {
            minlength = len;
            minlengthdefined = true;
          }
          if (maxlength < len)
          {
            maxlength = len;
          }
          sumlength += (uint64_t) len;
          numofseq++;
          if (arguments.dodistlen)
          {
            gt_disc_distri_add(distseqlen,len/arguments.bucketsize);
          }
          if (arguments.docstats)
          {
            gt_assembly_stats_calculator_add(asc, len);
          }
        }
        if (arguments.doastretch)
        {
          countA += accumulateastretch(distastretch,sequence,len);
        }
      }
      if (arguments.verbose)
      {
        gt_progressbar_stop();
      }
      gt_seq_iterator_delete(seqit);
    }
  }
  gt_str_array_delete(files);
  if (!had_err && arguments.dodistlen)
  {
    printf("# " Formatuint64_t " sequences of average length %.2f\n",
             PRINTuint64_tcast(numofseq),(double) sumlength/numofseq);
    printf("# total length " Formatuint64_t "\n",
             PRINTuint64_tcast(sumlength));
    printf("# minimum length %lu\n",minlength);
    printf("# maximum length %lu\n",maxlength);
    if (arguments.binarydistlen)
    {
      FILE *distlenoutfile;
      GtStr *distlenoutfilename = gt_str_new_cstr(argv[parsed_args]);
      gt_str_append_cstr(distlenoutfilename, GT_SEQSTAT_BINARY_DISTLEN_SUFFIX);
      distlenoutfile = gt_fa_fopen(gt_str_get(distlenoutfilename), "wb", err);
      if (distlenoutfile == NULL)
        had_err = -1;
      else
      {
        /*@ignore@*/
        gt_disc_distri_foreach(distseqlen, showdistseqlenbinary,
            distlenoutfile);
        /*@end@*/
        printf("# distribution of sequence length written to file: %s\n",
          gt_str_get(distlenoutfilename));
        gt_fa_fclose(distlenoutfile);
      }
      gt_str_delete(distlenoutfilename);
    }
    else
    {
      printf("# distribution of sequence length in buckets of size %u\n",
          arguments.bucketsize);
      gt_disc_distri_foreach(distseqlen, showdistseqlen,
          &(arguments.bucketsize));
    }
    gt_disc_distri_delete(distseqlen);
  }
  if (!had_err && arguments.docstats)
  {
    asc_logger = gt_logger_new(true, GT_LOGGER_DEFLT_PREFIX, stdout);
    gt_assembly_stats_calculator_show(asc, asc_logger);
    gt_logger_delete(asc_logger);
  }
  if (asc != NULL)
    gt_assembly_stats_calculator_delete(asc);
  if (!had_err && arguments.doastretch)
  {
    processastretches(distastretch,countA);
    gt_disc_distri_delete(distastretch);
  }
  return had_err;
}
