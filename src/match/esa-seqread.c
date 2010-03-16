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

Sequentialsuffixarrayreader *newSequentialsuffixarrayreaderfromfile(
                                  const GtStr *indexname,
                                  unsigned int demand,
                                  GT_UNUSED Sequentialaccesstype seqactype,
                                  GtError *err)
{
  Sequentialsuffixarrayreader *ssar;

  ALLOCASSIGNSPACE(ssar,NULL,Sequentialsuffixarrayreader,1);
  ALLOCASSIGNSPACE(ssar->suffixarray,NULL,Suffixarray,1);
  if (mapsuffixarray (ssar->suffixarray,
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
  ssar->nextlcptabindex = (Seqpos) 1;
  ssar->largelcpindex = 0;
  ssar->numberofsuffixes = totallength+1;
  return ssar;
}

void freeSequentialsuffixarrayreader(Sequentialsuffixarrayreader **ssar)
{
  if ((*ssar)->suffixarray != NULL)
  {
    freesuffixarray((*ssar)->suffixarray);
    FREESPACE((*ssar)->suffixarray);
  }
  FREESPACE(*ssar);
}

int nextSequentialsuftabvalue(Seqpos *currentsuffix,
                              Sequentialsuffixarrayreader *ssar,
                              GT_UNUSED GtError *err)
{
  *currentsuffix = ssar->suffixarray->suftab[ssar->nextsuftabindex++];
  return 1;
}

const GtEncodedsequence *encseqSequentialsuffixarrayreader(
                          const Sequentialsuffixarrayreader *ssar)
{
  return ssar->suffixarray->encseq;
}

Readmode readmodeSequentialsuffixarrayreader(
                          const Sequentialsuffixarrayreader *ssar)
{
  return ssar->suffixarray->readmode;
}

#else

 DECLAREREADFUNCTION(Seqpos);

 DECLAREREADFUNCTION(GtUchar);

 DECLAREREADFUNCTION(Largelcpvalue);

 struct Sequentialsuffixarrayreader
{
  Suffixarray *suffixarray;
  Seqpos numberofsuffixes,
         nextsuftabindex, /* for SEQ_mappedboth | SEQ_suftabfrommemory */
         nextlcptabindex, /* for SEQ_mappedboth */
         largelcpindex;   /* SEQ_mappedboth */
  Sequentialaccesstype seqactype;
  Lcpvalueiterator *lvi;
  const Seqpos *suftab;
  const GtEncodedsequence *encseq;
  Readmode readmode;
};

Sequentialsuffixarrayreader *newSequentialsuffixarrayreaderfromfile(
                                        const GtStr *indexname,
                                        unsigned int demand,
                                        Sequentialaccesstype seqactype,
                                        GtError *err)
{
  Sequentialsuffixarrayreader *ssar;

  ALLOCASSIGNSPACE(ssar,NULL,Sequentialsuffixarrayreader,1);
  ALLOCASSIGNSPACE(ssar->suffixarray,NULL,Suffixarray,1);
  gt_assert(seqactype == SEQ_mappedboth || seqactype == SEQ_scan);
  if (((seqactype == SEQ_mappedboth)
         ? mapsuffixarray : streamsuffixarray)(ssar->suffixarray,
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
  ssar->nextlcptabindex = (Seqpos) 1;
  ssar->largelcpindex = 0;
  ssar->seqactype = seqactype;
  ssar->suftab = NULL;
  ssar->encseq = ssar->suffixarray->encseq;
  ssar->readmode = ssar->suffixarray->readmode;
  ssar->numberofsuffixes = gt_encodedsequence_total_length(ssar->encseq) + 1;
  ssar->lvi = NULL;
  return ssar;
}

Sequentialsuffixarrayreader *newSequentialsuffixarrayreaderfromRAM(
                                        const GtEncodedsequence *encseq,
                                        Readmode readmode)
{
  Sequentialsuffixarrayreader *ssar;

  ALLOCASSIGNSPACE(ssar,NULL,Sequentialsuffixarrayreader,1);
  ssar->lvi = newLcpvalueiterator(encseq,readmode);
  ssar->suffixarray = NULL;
  ssar->nextlcptabindex = (Seqpos) 1; /* not required here */
  ssar->largelcpindex = 0; /* not required here */
  ssar->seqactype = SEQ_suftabfrommemory;
  ssar->readmode = readmode;
  ssar->encseq = encseq;
  return ssar;
}

void updateSequentialsuffixarrayreaderfromRAM(
                    Sequentialsuffixarrayreader *ssar,
                    const Seqpos *suftab,
                    bool firstpage,
                    Seqpos numberofsuffixes)
{
  ssar->nextsuftabindex = 0;
  ssar->suftab = suftab;
  ssar->numberofsuffixes = numberofsuffixes;
  if (firstpage)
  {
    (void) nextLcpvalueiterator(ssar->lvi,
                                true,
                                suftab,
                                numberofsuffixes);
  }
}

void freeSequentialsuffixarrayreader(Sequentialsuffixarrayreader **ssar)
{
  if ((*ssar)->suffixarray != NULL)
  {
    freesuffixarray((*ssar)->suffixarray);
    FREESPACE((*ssar)->suffixarray);
  }
  if ((*ssar)->lvi != NULL)
  {
    freeLcpvalueiterator(&(*ssar)->lvi);
  }
  FREESPACE(*ssar);
}

int nextSequentiallcpvalue(Seqpos *currentlcp,
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
    *currentlcp = nextLcpvalueiterator(ssar->lvi,
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
      *currentlcp = (Seqpos) tmpsmalllcpvalue;
    }
  }
  return 1;
}

int nextSequentialsuftabvalue(Seqpos *currentsuffix,
                              Sequentialsuffixarrayreader *ssar)
{
  if (ssar->seqactype == SEQ_scan)
  {
    return readnextSeqposfromstream(currentsuffix,
                                    &ssar->suffixarray->suftabstream);
  }
  if (ssar->seqactype == SEQ_mappedboth)
  {
    *currentsuffix = ssar->suffixarray->suftab[ssar->nextsuftabindex++];
    return 1;
  }
  *currentsuffix = ssar->suftab[ssar->nextsuftabindex++];
  return 1;
}

const GtEncodedsequence *encseqSequentialsuffixarrayreader(
                          const Sequentialsuffixarrayreader *ssar)
{
  return ssar->encseq;
}

Readmode readmodeSequentialsuffixarrayreader(
                          const Sequentialsuffixarrayreader *ssar)
{
  return ssar->readmode;
}
#endif /* ifdef INLINEDSequentialsuffixarrayreader */

const Seqpos *suftabSequentialsuffixarrayreader(
              const Sequentialsuffixarrayreader *ssar)
{
  gt_assert(ssar->seqactype != SEQ_scan);
  if (ssar->seqactype == SEQ_mappedboth)
  {
    return ssar->suffixarray->suftab;
  }
  return ssar->suftab;
}

const Suffixarray *suffixarraySequentialsuffixarrayreader(
              const Sequentialsuffixarrayreader *ssar)
{
  return ssar->suffixarray;
}
