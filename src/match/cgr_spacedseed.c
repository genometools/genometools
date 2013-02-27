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
#include "core/seq_iterator_sequence_buffer_api.h"
#include "cgr_spacedseed.h"
#include "sarr-def.h"
#include "core/intbits.h"
#include "idx-limdfs.h"
#include "absdfstrans-def.h"
#include "spaced-seeds.h"
#include "iter-window.h"
#include "esa-map.h"

typedef struct
{
  unsigned long seedwidth, numofonepositions;
  unsigned char *onepositions;
  GtBitsequence seedbitvector;
} Spacedseed;

static Spacedseed *spacedseed_new(const char *seedstring, GtError *err)
{
  Spacedseed *spse;
  unsigned long idx, fillindex;
  bool haserr = false;

  spse = gt_malloc(sizeof (*spse));
  spse->seedwidth = spse->numofonepositions = 0;
  spse->seedbitvector = 0;
  spse->onepositions = NULL;
  for (idx = 0; seedstring[idx] != '\0'; idx++)
  {
    if (spse->seedwidth == (unsigned long) GT_INTWORDSIZE)
    {
      gt_error_set(err,"space seed cannot be longer than %u",
                       (unsigned int) GT_INTWORDSIZE);
      haserr = true;
      break;
    }
    spse->seedwidth++;
    if (seedstring[idx] == '1')
    {
      spse->numofonepositions++;
      spse->seedbitvector |= GT_ITHBIT(idx);
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
                                        const GtUchar *qptr,
                                        const Spacedseed *spse)
{
  gt_indexbasedspacedseeds(limdfsresources,
                        qptr,
                        spse->seedbitvector,
                        spse->seedwidth,
                        dfst);
}

static void singlequerymatchspacedseed(Limdfsresources *limdfsresources,
                                       const AbstractDfstransformer *dfst,
                                       const GtUchar *query,
                                       unsigned long querylen,
                                       const Spacedseed *spse)
{
  const GtUchar *qptr;
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

static void cgr_showmatch(GT_UNUSED void *processinfo,
                      const GtIdxMatch *match)
{
  printf("%lu\t",match->dblen);
  printf("%lu\n",match->dbstartpos);
}

#ifdef WITHONLINE
static void onlinespacedseedsearch(const GtEncseq *encseq,
                                   const Spacedseed *spse,
                                   const GtUchar *qptr,qp)
{
  Windowiterator *wit;
  const GtUchar *buffer;
  unsigned long currentpos, totallength;
  unsigned long firstpos, windowschecked = 0;
  Bitsequence bitmask;
  bool matched;

  totallength = gt_encseq_total_length(encseq);
  wit = gt_windowiterator_new(encseq,spse->seedwidth,0,totallength);
  while (true)
  {
    buffer = gt_windowiterator_next(&currentpos,&firstpos,wit);
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
  gt_windowiterator_delete(wit);
}
#endif

int gt_matchspacedseed(bool withesa,
                    bool docompare,
                    const char *inputindex,
                    const GtStrArray *queryfilenames,
                    bool verbose,
                    GtError *err)
{
  Genericindex *genericindex = NULL;
  bool haserr = false;
  Spacedseed *spse;
  GtLogger *logger;

  logger = gt_logger_new(verbose, GT_LOGGER_DEFLT_PREFIX, stdout);
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
    genericindex = genericindex_new(inputindex,withesa,
                                    withesa && docompare,false,false,
                                    0,logger,err);
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
          0,gt_encseq_total_length(genericindex->suffixarray->encseq));
    iteroverallwords2(genericindex->suffixarray->encseq,spse->seedwidth,
            0,
            gt_encseq_total_length(genericindex->suffixarray->encseq));
    */
  }
  if (!haserr)
  {
    GtSeqIterator *seqit;
    const GtUchar *query;
    unsigned long querylen;
    char *desc = NULL;
    uint64_t unitnum;
    int retval;
    Limdfsresources *limdfsresources = NULL;
    const AbstractDfstransformer *dfst;
    const GtEncseq *encseq;

    dfst = gt_spse_AbstractDfstransformer();
    gt_assert(genericindex != NULL);
    limdfsresources = gt_newLimdfsresources(genericindex,
                                         true,
                                         0,
                                         (unsigned long) GT_INTWORDSIZE,
                                         false, /* keepexpandedonstack */
                                         cgr_showmatch,
                                         NULL, /* processmatch info */
                                         NULL, /* processresult */
                                         NULL, /* processresult info */
                                         dfst);
    encseq = genericindex_getencseq(genericindex);
    seqit = gt_seq_iterator_sequence_buffer_new(queryfilenames, err);
    if (!seqit)
      haserr = true;
    if (!haserr)
    {
      GtAlphabet *a = gt_encseq_alphabet(encseq);
      gt_seq_iterator_set_symbolmap(seqit, gt_alphabet_symbolmap(a));
      for (unitnum = 0; /* Nothing */; unitnum++)
      {
        retval = gt_seq_iterator_next(seqit,
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
      }
      if (limdfsresources != NULL)
      {
        gt_freeLimdfsresources(&limdfsresources,dfst);
      }
      gt_seq_iterator_delete(seqit);
    }
  }
  genericindex_delete(genericindex);
  spacedseed_delete(spse);
  gt_logger_delete(logger);
  return haserr ? -1 : 0;
}
