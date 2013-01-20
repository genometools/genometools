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
#include "esa-lcpval.h"
#include "lcpoverflow.h"
#include "esa-map.h"

Sequentialsuffixarrayreader *gt_newSequentialsuffixarrayreaderfromfile(
                                        const char *indexname,
                                        unsigned int demand,
                                        Sequentialaccesstype seqactype,
                                        GtLogger *logger,
                                        GtError *err)
{
  Sequentialsuffixarrayreader *ssar;

  ssar = gt_malloc(sizeof *ssar);
  ssar->suffixarray = gt_malloc(sizeof *ssar->suffixarray);
  gt_assert(seqactype == SEQ_mappedboth || seqactype == SEQ_scan);
  if ((seqactype == SEQ_mappedboth
         ? gt_mapsuffixarray : streamsuffixarray)(ssar->suffixarray,
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
  ssar->seqactype = seqactype;
  ssar->suftab = NULL;
  gt_assert(ssar->suffixarray != NULL);
  ssar->encseq = ssar->suffixarray->encseq;
  ssar->readmode = ssar->suffixarray->readmode;
  ssar->numberofsuffixes = gt_encseq_total_length(ssar->encseq) + 1;
  ssar->nonspecials = gt_encseq_total_length(ssar->encseq) -
                      gt_encseq_specialcharacters(ssar->encseq);
  ssar->lvi = NULL;
  return ssar;
}

Sequentialsuffixarrayreader *gt_newSequentialsuffixarrayreaderfromRAM(
                                        const GtEncseq *encseq,
                                        GtReadmode readmode)
{
  Sequentialsuffixarrayreader *ssar;

  ssar = gt_malloc(sizeof *ssar);
  ssar->lvi = gt_newLcpvalueiterator(encseq,readmode);
  ssar->suffixarray = NULL;
  ssar->nextlcptabindex = 1UL; /* not required here */
  ssar->largelcpindex = 0; /* not required here */
  ssar->seqactype = SEQ_suftabfrommemory;
  ssar->readmode = readmode;
  ssar->encseq = encseq;
  ssar->numberofsuffixes = gt_encseq_total_length(encseq) + 1;
  ssar->nonspecials = gt_encseq_total_length(encseq) -
                      gt_encseq_specialcharacters(encseq);
  return ssar;
}

void gt_updateSequentialsuffixarrayreaderfromRAM(
                    Sequentialsuffixarrayreader *ssar,
                    const ESASuffixptr *suftab,
                    bool firstpage,
                    unsigned long numberofsuffixes)
{
  ssar->nextsuftabindex = 0;
  ssar->suftab = suftab;
  ssar->numberofsuffixes = numberofsuffixes;
  if (firstpage)
  {
    (void) gt_nextLcpvalueiterator(ssar->lvi,
                                   true,
                                   suftab,
                                   numberofsuffixes);
  }
}

void gt_freeSequentialsuffixarrayreader(Sequentialsuffixarrayreader **ssar)
{
  if ((*ssar)->suffixarray != NULL)
  {
    gt_freesuffixarray((*ssar)->suffixarray);
    gt_free((*ssar)->suffixarray);
  }
  gt_freeLcpvalueiterator((*ssar)->lvi);
  gt_free(*ssar);
}

int gt_nextSequentiallcpvalue(unsigned long *currentlcp,
                              Sequentialsuffixarrayreader *ssar,
                              GtError *err)
{
  GtUchar tmpsmalllcpvalue;
  int retval;

  switch (ssar->seqactype)
  {
    case SEQ_scan:
      retval = gt_readnextfromstream_GtUchar(&tmpsmalllcpvalue,
                                         &ssar->suffixarray->lcptabstream);
      if (retval > 0)
      {
        if (tmpsmalllcpvalue < LCPOVERFLOW)
        {
          *currentlcp = (unsigned long) tmpsmalllcpvalue;
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
      break;
    case SEQ_mappedboth:
      if (ssar->nextlcptabindex < ssar->numberofsuffixes)
      {
        tmpsmalllcpvalue = ssar->suffixarray->lcptab[ssar->nextlcptabindex++];
        if (tmpsmalllcpvalue < LCPOVERFLOW)
        {
          *currentlcp = (unsigned long) tmpsmalllcpvalue;
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
      break;
    case SEQ_suftabfrommemory:
      if (ssar->nextlcptabindex < ssar->numberofsuffixes)
      {
        *currentlcp = gt_nextLcpvalueiterator(ssar->lvi,
                                              true,
                                              ssar->suftab,
                                              ssar->numberofsuffixes);
        ssar->nextlcptabindex++;
      } else
      {
        return 0;
      }
      break;
  }
  return 1;
}

int gt_nextSequentialsuftabvalue(unsigned long *currentsuffix,
                                 Sequentialsuffixarrayreader *ssar)
{
  if (ssar->seqactype == SEQ_scan)
  {
#ifdef _LP64
    if (ssar->suffixarray->suftabstream_GtUlong.fp != NULL)
    {
      return gt_readnextfromstream_GtUlong(currentsuffix,
                                 &ssar->suffixarray->suftabstream_GtUlong);
    } else
    {
      uint32_t readvalue = 0;
      int ret = gt_readnextfromstream_uint32_t(
                                 &readvalue,
                                 &ssar->suffixarray->suftabstream_uint32_t);
      *currentsuffix = (unsigned long) readvalue;
      return ret;
    }
#else
     return gt_readnextfromstream_GtUlong(currentsuffix,
                                      &ssar->suffixarray->suftabstream_GtUlong);
#endif
  }
  if (ssar->seqactype == SEQ_mappedboth)
  {
    *currentsuffix = ESASUFFIXPTRGET(ssar->suffixarray->suftab,
                                     ssar->nextsuftabindex++);
    return 1;
  }
  *currentsuffix = ESASUFFIXPTRGET(ssar->suftab,ssar->nextsuftabindex++);
  return 1;
}

const GtEncseq *gt_encseqSequentialsuffixarrayreader(
                          const Sequentialsuffixarrayreader *ssar)
{
  return ssar->encseq;
}

GtReadmode gt_readmodeSequentialsuffixarrayreader(
                          const Sequentialsuffixarrayreader *ssar)
{
  return ssar->readmode;
}

const ESASuffixptr *gt_suftabSequentialsuffixarrayreader(
                        const Sequentialsuffixarrayreader *ssar)
{
  gt_assert(ssar->seqactype != SEQ_scan);
  if (ssar->seqactype == SEQ_mappedboth)
  {
    return ssar->suffixarray->suftab;
  }
  return ssar->suftab;
}

const Suffixarray *gt_suffixarraySequentialsuffixarrayreader(
              const Sequentialsuffixarrayreader *ssar)
{
  return ssar->suffixarray;
}

unsigned long gt_Sequentialsuffixarrayreader_nonspecials(
                          const Sequentialsuffixarrayreader *ssar)
{
  return ssar->nonspecials;
}

unsigned long gt_Sequentialsuffixarrayreader_totallength(
              const Sequentialsuffixarrayreader *ssar)
{
  gt_assert(ssar->numberofsuffixes > 0);
  return ssar->numberofsuffixes - 1;
}

unsigned int gt_Sequentialsuffixarrayreader_prefixlength(
              const Sequentialsuffixarrayreader *ssar)
{
  return ssar->suffixarray->prefixlength;
}

GtBcktab *gt_Sequentialsuffixarrayreader_bcktab(
              const Sequentialsuffixarrayreader *ssar)
{
  return ssar->suffixarray->bcktab;
}
