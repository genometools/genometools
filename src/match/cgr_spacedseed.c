/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include "core/unused_api.h"
#include "core/ma_api.h"
#include "core/seqiterator.h"
#include "cgr_spacedseed.h"
#include "sarr-def.h"
#include "intbits.h"
#include "idx-limdfs.h"
#include "absdfstrans-def.h"
#include "spaced-seeds.h"
#include "iter-window.h"
#include "esa-map.h"

typedef struct
{
  unsigned long seedwidth, numofonepositions;
  unsigned char *onepositions;
  Bitstring seedbitvector;
} Spacedseed;

static Spacedseed *spacedseed_new(const char *seedstring,GtError *err)
{
  Spacedseed *spse;
  unsigned long idx, fillindex;
  bool haserr = false;

  spse = gt_malloc(sizeof(*spse));
  spse->seedwidth = spse->numofonepositions = 0;
  spse->seedbitvector = 0;
  spse->onepositions = NULL;
  for (idx = 0; seedstring[idx] != '\0'; idx++)
  {
    if (spse->seedwidth == (unsigned long) INTWORDSIZE)
    {
      gt_error_set(err,"space seed cannot be longer than %u",
                       (unsigned int) INTWORDSIZE);
      haserr = true;
      break;
    }
    spse->seedwidth++;
    if (seedstring[idx] == '1')
    {
      spse->numofonepositions++;
      spse->seedbitvector |= ITHBIT(idx);
    }
  }
  if (!haserr)
  {
    spse->onepositions = gt_malloc(sizeof (*spse->onepositions) *
                                   spse->numofonepositions);
    for (fillindex = 0, idx = 0; seedstring[idx] != '\0'; idx++)
    {
      if (seedstring[idx] == '1')
      {
        spse->onepositions[fillindex++] = (unsigned char) idx;
      }
    }
  }
  if (haserr)
  {
    gt_free(spse);
    return NULL;
  }
  return spse;
}

static void spacedseed_delete(Spacedseed *spse)
{
  gt_free(spse->onepositions);
  gt_free(spse);
}

static void singlewindowmatchspacedseed(Limdfsresources *limdfsresources,
                                        const AbstractDfstransformer *dfst,
                                        const Uchar *qptr,
                                        const Spacedseed *spse)
{
  indexbasedspacedseeds(limdfsresources,
                        qptr,
                        spse->seedbitvector,
                        spse->seedwidth,
                        dfst);
}

static void singlequerymatchspacedseed(Limdfsresources *limdfsresources,
                                       const AbstractDfstransformer *dfst,
                                       const Uchar *query,
                                       unsigned long querylen,
                                       const Spacedseed *spse)
{
  const Uchar *qptr;
  unsigned long offset, skipvalue;

  if (spse->seedwidth > querylen)
  {
    return;
  }
  qptr = query;
  offset = 0;
  while (qptr <= query + querylen - spse->seedwidth)
  {
    skipvalue = containsspecialbytestring(qptr,offset,spse->seedwidth);
    if (skipvalue == spse->seedwidth)
    {
      offset = spse->seedwidth-1;
      singlewindowmatchspacedseed(limdfsresources,dfst,qptr,spse);
      qptr++;
    } else
    {
      offset = 0;
      qptr += (skipvalue+1);
    }
  }
}

static void showmatch(GT_UNUSED void *processinfo,
                      Seqpos dbstartpos,
                      Seqpos dblen,
                      GT_UNUSED const Uchar *dbsubstring,
                      GT_UNUSED unsigned long pprefixlen,
                      GT_UNUSED unsigned long distance)
{
  printf(FormatSeqpos "\t",PRINTSeqposcast(dblen));
  printf(FormatSeqpos "\n",PRINTSeqposcast(dbstartpos));
}

#ifdef WITHONLINE
static void onlinespacedseedsearch(const Encodedsequence *encseq,
                                   const Spacedseed *spse,
                                   const Uchar *qptr,qp)
{
  Windowiterator *wit;
  const Uchar *buffer;
  Seqpos currentpos, totallength;
  unsigned long firstpos, windowschecked = 0;
  Bitstring bitmask;
  bool matched;

  totallength = getencseqtotallength(encseq);
  wit = windowiterator_new(encseq,spse->seedwidth,0,totallength);
  while (true)
  {
    buffer = windowiterator_next(&currentpos,&firstpos,wit);
    if (buffer != NULL)
    {
      bitmask = FIRSTBIT;
      matched = true;
      for (idx=0; idx < spse->seedwidth; idx++)
      {
        if ((spse->seedbitvector & bitmask) && qptr[idx] != buffer[idx])
        {
          matched = false;
          break;
        }
        bitmask >>= 1;
      }
      if (matched)
      {
      }
    } else
    {
      break;
    }
  }
  windowiterator_delete(wit);
}
#endif

int matchspacedseed(bool withesa,
                    bool docompare,
                    const GtStr *str_inputindex,
                    const GtStrArray *queryfilenames,
                    bool verbose,
                    GtError *err)
{
  Genericindex *genericindex = NULL;
  bool haserr = false;
  Spacedseed *spse;
  Verboseinfo *verboseinfo;

  verboseinfo = newverboseinfo(verbose);
  /*
  spse = spacedseed_new("11011011000011011",err);
  */
  spse = spacedseed_new("111001001001010111",err);
  if (spse == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    genericindex = genericindex_new(str_inputindex,withesa,
                                    withesa && docompare,0,verboseinfo,err);
    if (genericindex == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr && docompare)
  {
    gt_assert(genericindex != NULL);
    gt_assert(spse != NULL);
    /*
    iteroverallwords(genericindex->suffixarray->encseq,spse->seedwidth,
                     0,getencseqtotallength(genericindex->suffixarray->encseq));
    iteroverallwords2(genericindex->suffixarray->encseq,spse->seedwidth,
                      0,
                      getencseqtotallength(genericindex->suffixarray->encseq));
    */
  }
  if (!haserr)
  {
    GtSeqIterator *seqit;
    const Uchar *query;
    unsigned long querylen;
    char *desc = NULL;
    uint64_t unitnum;
    int retval;
    Limdfsresources *limdfsresources = NULL;
    const AbstractDfstransformer *dfst;
    const Encodedsequence *encseq;

    dfst = spse_AbstractDfstransformer();
    gt_assert(genericindex != NULL);
    limdfsresources = newLimdfsresources(genericindex,
                           true,
                           0,
                           (unsigned long) INTWORDSIZE,
                           showmatch,
                           NULL, /* processmatch info */
                           NULL, /* processresult */
                           NULL, /* processresult info */
                           dfst);
    encseq = genericindex_getencseq(genericindex);
    seqit = gt_seqiterator_new(queryfilenames,
                               getencseqAlphabetsymbolmap(encseq),
                               true);
    for (unitnum = 0; /* Nothing */; unitnum++)
    {
      retval = gt_seqiterator_next(seqit,
                                   &query,
                                   &querylen,
                                   &desc,
                                   err);
      if (retval < 0)
      {
        haserr = true;
        break;
      }
      if (retval == 0)
      {
        break;
      }
      singlequerymatchspacedseed(limdfsresources,
                                 dfst,
                                 query,
                                 querylen,
                                 spse);
      gt_free(desc);
    }
    if (limdfsresources != NULL)
    {
      freeLimdfsresources(&limdfsresources,dfst);
    }
    gt_seqiterator_delete(seqit);
  }
  genericindex_delete(genericindex);
  spacedseed_delete(spse);
  freeverboseinfo(&verboseinfo);
  return haserr ? -1 : 0;
}
