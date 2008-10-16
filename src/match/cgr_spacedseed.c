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

#include "esa-map.pr"

static int inputesaorpckindex(Suffixarray **suffixarray,
                              void **packedindex,
                              const GtStr *indexname,
                              bool withesa,
                              bool alwayswithencseq,
                              GtError *err)
{
  unsigned int demand;
  Seqpos totallength;
  bool haserr = false;

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
  *suffixarray = gt_malloc(sizeof(**suffixarray));
  if (mapsuffixarray(*suffixarray,
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
    if (withesa && (*suffixarray)->readmode != Forwardmode)
    {
      gt_error_set(err,"using option -esa you can only process index "
                       "in forward mode");
      haserr = true;
    } else
    {
      if (!withesa && (*suffixarray)->readmode != Reversemode)
      {
        gt_error_set(err,"with option -pck you can only process index "
                         "in reverse mode");
        haserr = true;
      }
    }
  }
  if (!haserr && !withesa)
  {
    *packedindex = loadvoidBWTSeqForSA(indexname,
                                       *suffixarray,
                                       totallength, true, err);
    if (*packedindex == NULL)
    {
      haserr = true;
    }
  }
  return haserr ? -1 : 0;
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

static void singlewindowmatchspacedseed(GT_UNUSED const void *genericindex,
                                        GT_UNUSED const Uchar *qptr,
                                        GT_UNUSED const Spacedseed *spse)
{
  return;
}

static void singlequerymatchspacedseed(const void *genericindex,
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

int matchspacedseed(bool withesa,
                    const GtStr *str_inputindex,
                    const GtStrArray *queryfilenames,
                    GtError *err)
{
  Suffixarray *suffixarray = NULL;
  void *packedindex = NULL;
  bool haserr = false;
  Spacedseed *spse;

  spse = spacedseed_new("11011011000011011",err);
  if (spse == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (inputesaorpckindex(&suffixarray,
                           &packedindex,
                           str_inputindex,
                           withesa,
                           false,
                           err) != 0)
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

    gt_assert(suffixarray != NULL);
    seqit = gt_seqiterator_new(queryfilenames,
                               getsymbolmapAlphabet(suffixarray->alpha),
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
      singlequerymatchspacedseed(withesa ? suffixarray : packedindex,
                                 query,
                                 querylen,
                                 spse);
      gt_free(desc);
    }
    gt_seqiterator_delete(seqit);
  }
  freesuffixarray(suffixarray);
  gt_free(suffixarray);
  if (packedindex != NULL)
  {
    deletevoidBWTSeq(packedindex);
  }
  spacedseed_delete(spse);
  return haserr ? -1 : 0;
}
