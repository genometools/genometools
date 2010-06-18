/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2010      Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2007-2010 Center for Bioinformatics, University of Hamburg

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

#include <sys/types.h>
#include <sys/stat.h>
#include <limits.h>
#include "core/fileutils_api.h"
#include "core/ma.h"
#include "core/option.h"
#include "core/seqiterator_sequence_buffer.h"
#include "core/seqiterator_fastq.h"
#include "core/unused_api.h"
#include "core/versionfunc.h"
#include "core/xposix.h"
#include "core/progressbar.h"
#include "core/disc_distri.h"
#include "core/format64.h"
#include "match/stamp.h"
#include "tools/gt_seqstat.h"

typedef struct
{
  bool verbose,
       dodistlen,
       doastretch,
       docstats,
       showestimsize;
  unsigned int bucketsize;
} GtSeqstatArguments;

static GtOPrval parse_options(GtSeqstatArguments *arguments,
                              int *parsed_args,int argc,
                              const char **argv, GtError *err)
{
  GtOptionParser *op;
  GtOption *optionverbose, *optiondistlen, *optionbucketsize,
           *optioncstats, *optionastretch, *optionestimsize;
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
                                &arguments->bucketsize,100, 1);
  gt_option_imply(optionbucketsize, optiondistlen);
  gt_option_parser_add_option(op, optionbucketsize);

  optioncstats = gt_option_new_bool("contigs",
                                   "summary of contigs set statistics",
                                   &arguments->docstats,false);
  gt_option_parser_add_option(op, optioncstats);

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

#define NOF_NSTATS 2 /* N50, N80 */

typedef struct
{
  unsigned long       nvalue[NOF_NSTATS];
  unsigned long long  min[NOF_NSTATS];
  bool                done[NOF_NSTATS];
  char                *name[NOF_NSTATS];
  unsigned long long  current;
} Nstats;

static void calcNstats(unsigned long key, unsigned long long value,
                        void *data)
{
  Nstats *nstats = data;
  int i;
  nstats->current += (key*value);
  for (i = 0; i < NOF_NSTATS; i++)
  {
    if (!nstats->done[i])
    {
      if (nstats->current >= nstats->min[i])
      {
        nstats->done[i] = true;
        nstats->nvalue[i] = key;
      }
    }
  }
}

#define initNstat(INDEX, NAME, PORTION)\
    nstats.name[INDEX] = (NAME);\
    nstats.min[INDEX] = (sumlength * (PORTION) / 100);\
    nstats.nvalue[INDEX] = 0;\
    nstats.done[INDEX] = false

static void showcstats(uint64_t numofseq, uint64_t sumlength,
                        unsigned long minlength, unsigned long maxlength,
                        GtDiscDistri *distctglen)
{
  Nstats nstats;
  int i;

  initNstat(0, "N50", 50);
  initNstat(1, "N80", 80);
  nstats.current = 0;

  gt_disc_distri_foreach_in_reverse_order(distctglen, calcNstats,
                                          &nstats);
  printf("# number of contigs: "Formatuint64_t"\n",
         PRINTuint64_tcast(numofseq));
  printf("# total length:      "Formatuint64_t"\n",
         PRINTuint64_tcast(sumlength));
  printf("# average size:      %.2f\n", (double) sumlength/numofseq);
  printf("# longest contig:    %lu\n", maxlength);
  for (i = 0; i < NOF_NSTATS ; i++)
    printf("# %s:               %lu\n", nstats.name[i], nstats.nvalue[i]);
  printf("# smallest contig:   %lu\n", minlength);
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
  GtDiscDistri *distctglen = NULL;
  GtDiscDistri *distastretch = NULL;
  uint64_t numofseq = 0, sumlength = 0;
  unsigned long minlength = 0, maxlength = 0;
  unsigned long long countA = 0;
  bool minlengthdefined = false;
  GtSeqstatArguments arguments;

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
    seqit = gt_seqiterator_sequence_buffer_new(files, err);
    if (!seqit)
      had_err = -1;
    if (!had_err)
    {
      if (arguments.dodistlen)
      {
        distseqlen = gt_disc_distri_new();
      }
      if (arguments.docstats)
      {
        distctglen = gt_disc_distri_new();
      }
      if (arguments.doastretch)
      {
        distastretch = gt_disc_distri_new();
      }
      if (arguments.verbose)
      {
        gt_progressbar_start(gt_seqiterator_getcurrentcounter(seqit,
                                                            (unsigned long long)
                                                            totalsize),
                             (unsigned long long) totalsize);
      }
      while (true)
      {
        desc = NULL;
        had_err = gt_seqiterator_next(seqit, &sequence, &len, &desc, err);
        gt_free(desc);
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
            gt_disc_distri_add(distctglen,len);
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
      gt_seqiterator_delete(seqit);
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
    printf("# distribution of sequence length in buckets of size %u\n",
           arguments.bucketsize);
    gt_disc_distri_foreach(distseqlen, showdistseqlen,
                           &(arguments.bucketsize));
    gt_disc_distri_delete(distseqlen);
  }
  if (!had_err && arguments.docstats)
  {
    showcstats(numofseq, sumlength, minlength,
               maxlength, distctglen);
  }
  if (distctglen) gt_disc_distri_delete(distctglen);
  if (!had_err && arguments.doastretch)
  {
    processastretches(distastretch,countA);
    gt_disc_distri_delete(distastretch);
  }
  return had_err;
}
