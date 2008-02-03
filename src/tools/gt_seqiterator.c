/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <sys/types.h>
#include <sys/stat.h>
#include "libgtcore/fileutils.h"
#include "libgtcore/ma.h"
#include "libgtcore/option.h"
#include "libgtcore/seqiterator.h"
#include "libgtcore/versionfunc.h"
#include "libgtcore/xposix.h"
#include "libgtcore/progressbar.h"
#include "libgtcore/discdistri.h"
#include "libgtmatch/format64.h"
#include "tools/gt_seqiterator.h"

#define BUCKETSIZE 100

static OPrval parse_options(bool *verbose,bool *dodistlen,int *parsed_args,
                            int argc, const char **argv, Error *err)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;

  error_check(err);
  op = option_parser_new("[options] file [...]",
                         "Parse the supplied Fasta files.");
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");
  option= option_new_bool("v","be verbose",verbose,false);
  option_parser_add_option(op, option);
  option= option_new_bool("distlen","show distribution of sequence length",
                           dodistlen,false);
  option_parser_add_option(op, option);
  option_parser_set_min_args(op, 1U);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);
  option_parser_delete(op);
  return oprval;
}

static void showdistseqlen(unsigned long key, unsigned long long value,
                           /*@unused@*/ void *data)
{
  unsigned long distvalue;

  distvalue = (unsigned long) value; /* XXX: is this cast always OK? */
  printf("%lu--%lu %lu\n",
         BUCKETSIZE * key,
         BUCKETSIZE * (key+1) - 1,
         distvalue);
}

int gt_seqiterator(int argc, const char **argv, Error *err)
{
  StrArray *files;
  SeqIterator *seqit;
  const Uchar *sequence;
  char *desc;
  unsigned long len;
  int i, parsed_args, had_err;
  off_t totalsize;
  DiscDistri *distseqlen = NULL;
  bool verbose = false, dodistlen = false;
  uint64_t numofseq = 0, sumlength = 0;
  unsigned long minlength = 0, maxlength = 0;
  bool minlengthdefined = false;

  error_check(err);

  /* option parsing */
  switch (parse_options(&verbose,&dodistlen,&parsed_args, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  files = strarray_new();
  for (i = parsed_args; i < argc; i++)
  {
    strarray_add_cstr(files, argv[i]);
  }

  totalsize = files_estimate_total_size(files);
  printf("# estimated total size is " Formatuint64_t "\n",
            PRINTuint64_tcast(totalsize));
  seqit = seqiterator_new(files, NULL, true);
  if (dodistlen)
  {
    distseqlen = discdistri_new();
  }
  if (verbose)
  {
    progressbar_start(seqiterator_getcurrentcounter(seqit, (unsigned long long)
                                                           totalsize),
                                                           (unsigned long long)
                                                           totalsize);
  }
  while (true)
  {
    had_err = seqiterator_next(seqit, &sequence, &len, &desc, err);
    if (dodistlen)
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
      discdistri_add(distseqlen,len/BUCKETSIZE);
    }
    if (had_err != 1)
    {
      break;
    }
    ma_free(desc);
  }
  if (verbose)
  {
    progressbar_stop();
  }
  seqiterator_delete(seqit);
  strarray_delete(files);
  if (dodistlen)
  {
    printf("# " Formatuint64_t " sequences of average length %.2f\n",
             PRINTuint64_tcast(numofseq),(double) sumlength/numofseq);
    printf("# minimum length %lu\n",minlength);
    printf("# maximum length %lu\n",maxlength);
    printf("# distribution of sequence length in buckets of size %d\n",
           BUCKETSIZE);
    discdistri_foreach(distseqlen,showdistseqlen,NULL);
    discdistri_delete(distseqlen);
  }
  return had_err;
}
