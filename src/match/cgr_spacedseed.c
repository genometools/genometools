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
#include "cgr_spacedseed.h"
#include "sarr-def.h"
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

int matchspacedseed(bool withesa,
                    const GtStr *str_inputindex,
                    GT_UNUSED const GtStrArray *queryfilenames,
                    GtError *err)
{
  Suffixarray *suffixarray = NULL;
  void *packedindex = NULL;
  bool haserr = false;

  if (inputesaorpckindex(&suffixarray,
                         &packedindex, 
                         str_inputindex,
                         withesa,
                         false,
                         err) != 0)
  {
    haserr = true;
  }
  freesuffixarray(suffixarray);
  gt_free(suffixarray);
  if (packedindex != NULL)
  {
    deletevoidBWTSeq(packedindex);
  }
  return haserr ? -1 : 0;
}
