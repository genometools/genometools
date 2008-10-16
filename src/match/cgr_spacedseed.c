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
#include "alphadef.h"
#include "sarr-def.h"
#include "intbits.h"
#include "eis-voiditf.h"
#include "idx-limdfs.h"
#include "absdfstrans-def.h"
#include "spaced-seeds.h"

#include "esa-map.pr"

typedef struct
{
  Suffixarray *suffixarray;
  void *packedindex;
  bool withesa;
  const Matchbound **mbtab; /* only relevant for packedindex */
  unsigned int maxdepth;    /* maximaldepth of boundaries */
} Genericindex;

static void genericindex_delete(Genericindex *genericindex)
{
  if (genericindex == NULL)
  {
    return;
  }
  freesuffixarray(genericindex->suffixarray);
  gt_free(genericindex->suffixarray);
  if (genericindex->packedindex != NULL)
  {
    deletevoidBWTSeq(genericindex->packedindex);
  }
  gt_free(genericindex);
}

static Genericindex *genericindex_new(const GtStr *indexname,
                                      bool withesa,
                                      bool alwayswithencseq,
                                      GtError *err)
{
  unsigned int demand;
  Seqpos totallength;
  bool haserr = false;
  Genericindex *genericindex;

  genericindex = gt_malloc(sizeof(*genericindex));
  if (withesa)
  {
    demand = SARR_ESQTAB | SARR_SUFTAB;
  } else
  {
    if (alwayswithencseq)
    {
      demand = SARR_ESQTAB;
    } else
    {
      demand = 0;
    }
  }
  genericindex->withesa = withesa;
  genericindex->suffixarray = gt_malloc(sizeof(*genericindex->suffixarray));
  if (mapsuffixarray(genericindex->suffixarray,
                     &totallength,
                     demand,
                     indexname,
                     NULL,
                     err) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (withesa && genericindex->suffixarray->readmode != Forwardmode)
    {
      gt_error_set(err,"using option -esa you can only process index "
                       "in forward mode");
      haserr = true;
    } else
    {
      if (!withesa && genericindex->suffixarray->readmode != Reversemode)
      {
        gt_error_set(err,"with option -pck you can only process index "
                         "in reverse mode");
        haserr = true;
      }
    }
  }
  genericindex->packedindex = NULL;
  genericindex->mbtab = NULL;
  genericindex->maxdepth = 0;
  if (!haserr && !withesa)
  {
    genericindex->packedindex = loadvoidBWTSeqForSA(indexname,
                                                    genericindex->suffixarray,
                                                    totallength, true, err);
    if (genericindex->packedindex == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr && !withesa)
  {
    genericindex->mbtab = bwtseq2mbtab(genericindex->packedindex);
    genericindex->maxdepth = bwtseq2maxdepth(genericindex->packedindex);
  }
  if (haserr)
  {
    genericindex_delete(genericindex);
    return NULL;
  }
  return genericindex;
}

typedef struct
{
  unsigned long seedweight, numofonepositions;
  unsigned char *onepositions;
  Bitstring seedbitvector;
} Spacedseed;

static Spacedseed *spacedseed_new(const char *seedstring,GtError *err)
{
  Spacedseed *spse;
  unsigned long idx, fillindex;
  bool haserr = false;

  spse = gt_malloc(sizeof(*spse));
  spse->seedweight = spse->numofonepositions = 0;
  spse->seedbitvector = 0;
  spse->onepositions = NULL;
  for (idx = 0; seedstring[idx] != '\0'; idx++)
  {
    if (spse->seedweight == (unsigned long) INTWORDSIZE)
    {
      gt_error_set(err,"space seed cannot be longer than %u",
                       (unsigned int) INTWORDSIZE);
      haserr = true;
      break;
    }
    spse->seedweight++;
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

static void singlewindowmatchspacedseed(
                        GT_UNUSED const Genericindex *genericindex,
                        GT_UNUSED const Uchar *qptr,
                        GT_UNUSED const Spacedseed *spse)
{
  return;
}

static void singlequerymatchspacedseed(const Genericindex *genericindex,
                                       const Uchar *query,
                                       unsigned long querylen,
                                       const Spacedseed *spse)
{
  const Uchar *qptr;
  unsigned long offset, skipvalue;

  if (spse->seedweight > querylen)
  {
    return;
  }
  qptr = query;
  offset = 0;
  while (qptr <= query + querylen - spse->seedweight)
  {
    skipvalue = containsspecialbytestring(qptr,offset,spse->seedweight);
    if (skipvalue == spse->seedweight)
    {
      offset = spse->seedweight-1;
      singlewindowmatchspacedseed(genericindex,qptr,spse);
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
                      GT_UNUSED unsigned long pprefixlen)
{
  printf(FormatSeqpos "\t",PRINTSeqposcast(dblen));
  printf(FormatSeqpos "\n",PRINTSeqposcast(dbstartpos));
}

int matchspacedseed(bool withesa,
                    const GtStr *str_inputindex,
                    const GtStrArray *queryfilenames,
                    GtError *err)
{
  Genericindex *genericindex = NULL;
  bool haserr = false;
  Spacedseed *spse;

  spse = spacedseed_new("11011011000011011",err);
  if (spse == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    genericindex = genericindex_new(str_inputindex,withesa,false,err);
    if (genericindex == NULL)
    {
      haserr = true;
    }
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

    dfst = spse_AbstractDfstransformer();
    gt_assert(genericindex != NULL);
    gt_assert(genericindex->suffixarray != NULL);
    limdfsresources = newLimdfsresources(
                           withesa ? genericindex->suffixarray
                                   : genericindex->packedindex,
                           genericindex->mbtab,
                           genericindex->maxdepth,
                           genericindex->suffixarray->encseq,
                           withesa,
                           true,
                           0,
                           getmapsizeAlphabet(genericindex->suffixarray->alpha),
                           getencseqtotallength(genericindex->suffixarray
                                                            ->encseq),
                           (unsigned long) INTWORDSIZE,
                           showmatch,
                           NULL, /* processmatch info */
                           NULL, /* processresult */
                           NULL, /* processresult info */
                           dfst);
    seqit = gt_seqiterator_new(queryfilenames,
                               getsymbolmapAlphabet(genericindex->suffixarray
                                                                ->alpha),
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
      singlequerymatchspacedseed(genericindex,
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
  return haserr ? -1 : 0;
}
