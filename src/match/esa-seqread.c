/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include <limits.h>
#include "core/unused_api.h"
#include "core/ma_api.h"
#include "sarr-def.h"
#include "esa-seqread.h"
#include "lcpoverflow.h"
#include "esa-map.h"

Sequentialsuffixarrayreader *gt_newSequentialsuffixarrayreaderfromfile(
                                        const char *indexname,
                                        unsigned int demand,
                                        bool scanfile,
                                        GtLogger *logger,
                                        GtError *err)
{
  Sequentialsuffixarrayreader *ssar;

  ssar = gt_malloc(sizeof *ssar);
  ssar->suffixarray = gt_malloc(sizeof *ssar->suffixarray);
  if ((scanfile ? streamsuffixarray : gt_mapsuffixarray)(ssar->suffixarray,
                                                         demand,
                                                         indexname,
                                                         logger,
                                                         err) != 0)
  {
    gt_free(ssar->suffixarray);
    gt_free(ssar);
    return NULL;
  }
  ssar->nextsuftabindex = 0;
  ssar->nextlcptabindex = 1UL;
  ssar->largelcpindex = 0;
  ssar->scanfile = scanfile;
  ssar->suftab = NULL;
  gt_assert(ssar->suffixarray != NULL);
  ssar->encseq = ssar->suffixarray->encseq;
  ssar->readmode = ssar->suffixarray->readmode;
  ssar->numberofsuffixes = gt_encseq_total_length(ssar->encseq) + 1;
  ssar->nonspecials = gt_encseq_total_length(ssar->encseq) -
                      gt_encseq_specialcharacters(ssar->encseq);
  ssar->extrainfo = NULL;
  return ssar;
}

void gt_freeSequentialsuffixarrayreader(Sequentialsuffixarrayreader **ssar)
{
  if ((*ssar)->suffixarray != NULL)
  {
    gt_freesuffixarray((*ssar)->suffixarray);
    gt_free((*ssar)->suffixarray);
  }
  gt_free(*ssar);
}

int gt_nextSequentiallcpvalue(GtUword *currentlcp,
                              Sequentialsuffixarrayreader *ssar,
                              GtError *err)
{
  GtUchar tmpsmalllcpvalue;
  int retval;

  gt_assert(ssar != NULL);
  if (ssar->scanfile)
  {
    retval = gt_readnextfromstream_GtUchar(&tmpsmalllcpvalue,
                                           &ssar->suffixarray->lcptabstream);
    if (retval > 0)
    {
      if (tmpsmalllcpvalue < (GtUchar) LCPOVERFLOW)
      {
        *currentlcp = (GtUword) tmpsmalllcpvalue;
      } else
      {
        Largelcpvalue tmpexception;

        retval = gt_readnextfromstream_Largelcpvalue(&tmpexception,
                                          &ssar->suffixarray->llvtabstream);
        if (retval == 0)
        {
          gt_error_set(err,"file %s: line %d: unexpected end of file when "
                           "reading llvtab",__FILE__,__LINE__);
          return -1;
        }
        *currentlcp = tmpexception.value;
      }
    } else
    {
      return 0;
    }
  } else
  {
    if (ssar->nextlcptabindex < ssar->numberofsuffixes)
    {
      tmpsmalllcpvalue = ssar->suffixarray->lcptab[ssar->nextlcptabindex++];
      if (tmpsmalllcpvalue < (GtUchar) LCPOVERFLOW)
      {
        *currentlcp = (GtUword) tmpsmalllcpvalue;
      } else
      {
        gt_assert(ssar->suffixarray->llvtab[ssar->largelcpindex].position ==
                  ssar->nextlcptabindex-1);
        *currentlcp = ssar->suffixarray->llvtab[ssar->largelcpindex++].value;
      }
    } else
    {
      return 0;
    }
  }
  return 1;
}

int gt_nextSequentialsuftabvalue(GtUword *currentsuffix,
                                 Sequentialsuffixarrayreader *ssar)
{
  gt_assert(ssar != NULL);
  if (ssar->scanfile)
  {
#if defined (_LP64) || defined (_WIN64)
    if (ssar->suffixarray->suftabstream_GtUword.fp != NULL)
    {
      return gt_readnextfromstream_GtUword(currentsuffix,
                                 &ssar->suffixarray->suftabstream_GtUword);
    } else
    {
      uint32_t readvalue = 0;
      int ret = gt_readnextfromstream_uint32_t(
                                 &readvalue,
                                 &ssar->suffixarray->suftabstream_uint32_t);
      *currentsuffix = (GtUword) readvalue;
      return ret;
    }
#else
     return gt_readnextfromstream_GtUword(currentsuffix,
                                      &ssar->suffixarray->suftabstream_GtUword);
#endif
  }
  if (ssar->scanfile)
  {
    *currentsuffix = ESASUFFIXPTRGET(ssar->suftab,ssar->nextsuftabindex++);
  } else
  {
    *currentsuffix = ESASUFFIXPTRGET(ssar->suffixarray->suftab,
                                     ssar->nextsuftabindex++);
  }
  return 1;
}

const GtEncseq *gt_encseqSequentialsuffixarrayreader(
                          const Sequentialsuffixarrayreader *ssar)
{
  gt_assert(ssar != NULL);
  return ssar->encseq;
}

GtReadmode gt_readmodeSequentialsuffixarrayreader(
                          const Sequentialsuffixarrayreader *ssar)
{
  gt_assert(ssar != NULL);
  return ssar->readmode;
}

const ESASuffixptr *gt_suftabSequentialsuffixarrayreader(
                        const Sequentialsuffixarrayreader *ssar)
{
  gt_assert(ssar != NULL);
  if (ssar->scanfile)
  {
    return ssar->suftab;
  }
  return ssar->suffixarray->suftab;
}

const Suffixarray *gt_suffixarraySequentialsuffixarrayreader(
              const Sequentialsuffixarrayreader *ssar)
{
  gt_assert(ssar != NULL);
  return ssar->suffixarray;
}

GtUword gt_Sequentialsuffixarrayreader_nonspecials(
                          const Sequentialsuffixarrayreader *ssar)
{
  gt_assert(ssar != NULL);
  return ssar->nonspecials;
}

GtUword gt_Sequentialsuffixarrayreader_totallength(
              const Sequentialsuffixarrayreader *ssar)
{
  gt_assert(ssar != NULL && ssar->numberofsuffixes > 0);
  return ssar->numberofsuffixes - 1;
}

GtUword gt_Sequentialsuffixarrayreader_maxbranchdepth(
              const Sequentialsuffixarrayreader *ssar)
{
  gt_assert(ssar != NULL && ssar->suffixarray != NULL &&
            ssar->suffixarray->maxbranchdepth.defined);
  return ssar->suffixarray->maxbranchdepth.valueunsignedlong;
}

unsigned int gt_Sequentialsuffixarrayreader_prefixlength(
              const Sequentialsuffixarrayreader *ssar)
{
  gt_assert(ssar != NULL && ssar->suffixarray != NULL);
  return ssar->suffixarray->prefixlength;
}
