/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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
#include "alphadef.h"
#include "sarr-def.h"
#include "intbits.h"
#include "eis-voiditf.h"
#include "idx-limdfs.h"
#include "idxlocali.h"
#include "idxlocalidp.h"
#include "absdfstrans-def.h"
#include "esa-map.h"

typedef struct
{
  Suffixarray *suffixarray;
  Seqpos totallength;
  void *packedindex;
  bool withesa;
  const Mbtab **mbtab;      /* only relevant for packedindex */
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
                     demand,
                     indexname,
                     NULL,
                     err) != 0)
  {
    haserr = true;
    genericindex->totallength = 0;
  } else
  {
    genericindex->totallength = getencseqtotallength(genericindex->suffixarray
                                                                 ->encseq);
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
                                                    genericindex->totallength,
                                                    true, err);
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

int runidxlocali(const IdxlocaliOptions *arguments,GtError *err)
{
  Genericindex *genericindex = NULL;
  bool haserr = false;

  genericindex = genericindex_new(arguments->indexname,arguments->withesa,
                                  false,err);
  if (genericindex == NULL)
  {
    haserr = true;
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

    dfst = locali_AbstractDfstransformer();
    gt_assert(genericindex != NULL);
    gt_assert(genericindex->suffixarray != NULL);
    limdfsresources = newLimdfsresources(
                           arguments->withesa ? genericindex->suffixarray
                                              : genericindex->packedindex,
                           genericindex->mbtab,
                           genericindex->maxdepth,
                           genericindex->suffixarray->encseq,
                           arguments->withesa,
                           true,
                           0,
                           genericindex->totallength,
                           0, /* XXX maxpathlength? check */
                           showmatch,
                           NULL, /* processmatch info */
                           NULL, /* processresult */
                           NULL, /* processresult info */
                           dfst);
    seqit = gt_seqiterator_new(arguments->queryfiles,
                               getencseqAlphabetsymbolmap(genericindex->
                                                          suffixarray->encseq),
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
      indexbasedlocali(limdfsresources,
                       arguments->matchscore,
                       arguments->mismatchscore,
                       arguments->gapstart,
                       arguments->gapextend,
                       arguments->threshold,
                       query,
                       querylen,
                       dfst);
      gt_free(desc);
    }
    if (limdfsresources != NULL)
    {
      freeLimdfsresources(&limdfsresources,dfst);
    }
    gt_seqiterator_delete(seqit);
  }
  genericindex_delete(genericindex);
  return haserr ? -1 : 0;
}
