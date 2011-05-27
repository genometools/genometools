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
#include "sarr-def.h"
#include "spacedef.h"
#include "esa-seqread.h"
#include "esa-lcpval.h"
#include "lcpoverflow.h"
#include "esa-map.h"

#ifdef INLINEDSequentialsuffixarrayreader

Sequentialsuffixarrayreader *gt_newSequentialsuffixarrayreaderfromfile(
                                  const char *indexname,
                                  unsigned int demand,
                                  GT_UNUSED Sequentialaccesstype seqactype,
                                  GtError *err)
{
  Sequentialsuffixarrayreader *ssar;

  /* printf("INLINEDSequentialsuffixarrayreader=true\n"); */
  ALLOCASSIGNSPACE(ssar,NULL,Sequentialsuffixarrayreader,1);
  ALLOCASSIGNSPACE(ssar->suffixarray,NULL,Suffixarray,1);
  if (gt_mapsuffixarray (ssar->suffixarray,
                         demand,
                         indexname,
                         NULL,
                         err) != 0)
  {
    FREESPACE(ssar->suffixarray);
    FREESPACE(ssar);
    return NULL;
  }
  ssar->nextsuftabindex = 0;
  ssar->nextlcptabindex = 1UL;
  ssar->largelcpindex = 0;
  gt_assert(ssar->suffixarray != NULL);
  ssar->numberofsuffixes
    = gt_encseq_total_length(ssar->suffixarray->encseq) + 1;
  ssar->nonspecials = gt_encseq_total_length(ssar->suffixarray->encseq) -
                      gt_encseq_specialcharacters(ssar->suffixarray->encseq);
  return ssar;
}

void gt_freeSequentialsuffixarrayreader(Sequentialsuffixarrayreader **ssar)
{
  if ((*ssar)->suffixarray != NULL)
  {
    gt_freesuffixarray((*ssar)->suffixarray);
    FREESPACE((*ssar)->suffixarray);
  }
  FREESPACE(*ssar);
}

int gt_nextSequentialsuftabvalue(unsigned long *currentsuffix,
                              Sequentialsuffixarrayreader *ssar,
                              GT_UNUSED GtError *err)
{
  gt_assert(ssar != NULL && ssar->suffixarray != NULL &&
            ssar->suffixarray->suftab != NULL);
  *currentsuffix = ssar->suffixarray->suftab[ssar->nextsuftabindex++];
  return 1;
}

const GtEncseq *gt_encseqSequentialsuffixarrayreader(
                          const Sequentialsuffixarrayreader *ssar)
{
  gt_assert(ssar != NULL && ssar->suffixarray != NULL);
  return ssar->suffixarray->encseq;
}

GtReadmode gt_readmodeSequentialsuffixarrayreader(
                          const Sequentialsuffixarrayreader *ssar)
{
  gt_assert(ssar != NULL && ssar->suffixarray != NULL);
  return ssar->suffixarray->readmode;
}

const ESASuffixptr *gt_suftabSequentialsuffixarrayreader(
                          const Sequentialsuffixarrayreader *ssar)
{
  gt_assert(ssar != NULL && ssar->suffixarray != NULL);
  return ssar->suffixarray->suftab;
}

#else

struct Sequentialsuffixarrayreader
{
  Suffixarray *suffixarray;
  unsigned long nonspecials,
         numberofsuffixes,
         nextsuftabindex, /* for SEQ_mappedboth | SEQ_suftabfrommemory */
         nextlcptabindex, /* for SEQ_mappedboth */
         largelcpindex;   /* SEQ_mappedboth */
  Sequentialaccesstype seqactype;
  Lcpvalueiterator *lvi;
  const ESASuffixptr *suftab;
  const GtEncseq *encseq;
  GtReadmode readmode;
};

Sequentialsuffixarrayreader *gt_newSequentialsuffixarrayreaderfromfile(
                                        const char *indexname,
                                        unsigned int demand,
                                        Sequentialaccesstype seqactype,
                                        GtError *err)
{
  Sequentialsuffixarrayreader *ssar;

  ALLOCASSIGNSPACE(ssar,NULL,Sequentialsuffixarrayreader,1);
  ALLOCASSIGNSPACE(ssar->suffixarray,NULL,Suffixarray,1);
  gt_assert(seqactype == SEQ_mappedboth || seqactype == SEQ_scan);
  if (((seqactype == SEQ_mappedboth)
         ? gt_mapsuffixarray : streamsuffixarray)(ssar->suffixarray,
                                               demand,
                                               indexname,
                                               NULL,
                                               err) != 0)
  {
    FREESPACE(ssar->suffixarray);
    FREESPACE(ssar);
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

  ALLOCASSIGNSPACE(ssar,NULL,Sequentialsuffixarrayreader,1);
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
    FREESPACE((*ssar)->suffixarray);
  }
  if ((*ssar)->lvi != NULL)
  {
    gt_freeLcpvalueiterator(&(*ssar)->lvi);
  }
  FREESPACE(*ssar);
}

int gt_nextSequentiallcpvalue(unsigned long *currentlcp,
                           Sequentialsuffixarrayreader *ssar,
                           GtError *err)
{
  GtUchar tmpsmalllcpvalue;
  int retval;

  if (ssar->seqactype == SEQ_suftabfrommemory)
  {
    if (ssar->nextlcptabindex >= ssar->numberofsuffixes)
    {
      return 0;
    }
    *currentlcp = gt_nextLcpvalueiterator(ssar->lvi,
                                       true,
                                       ssar->suftab,
                                       ssar->numberofsuffixes);
    ssar->nextlcptabindex++;
  } else
  {
    if (ssar->seqactype == SEQ_mappedboth)
    {
      if (ssar->nextlcptabindex >= ssar->numberofsuffixes)
      {
        return 0;
      }
      tmpsmalllcpvalue = ssar->suffixarray->lcptab[ssar->nextlcptabindex++];
    } else
    {
      retval = readnextGtUcharfromstream(&tmpsmalllcpvalue,
                                       &ssar->suffixarray->lcptabstream);
      if (retval == 0)
      {
        return 0;
      }
    }
    if (tmpsmalllcpvalue == LCPOVERFLOW)
    {
      Largelcpvalue tmpexception;

      if (ssar->seqactype == SEQ_mappedboth)
      {
        gt_assert(ssar->suffixarray->llvtab[ssar->largelcpindex].position ==
               ssar->nextlcptabindex-1);
        *currentlcp = ssar->suffixarray->llvtab[ssar->largelcpindex++].value;
      } else
      {
        retval = readnextLargelcpvaluefromstream(
                                          &tmpexception,
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
      *currentlcp = (unsigned long) tmpsmalllcpvalue;
    }
  }
  return 1;
}

int gt_nextSequentialsuftabvalue(unsigned long *currentsuffix,
                              Sequentialsuffixarrayreader *ssar)
{
  if (ssar->seqactype == SEQ_scan)
  {
    return readnextGtUlongfromstream(currentsuffix,
                                    &ssar->suffixarray->suftabstream);
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
#endif /* ifdef INLINEDSequentialsuffixarrayreader */

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
