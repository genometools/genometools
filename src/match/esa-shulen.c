/*
  Copyright (c) 2010 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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
#include "core/array2dim_api.h"
#include "core/logger.h"
#include "esa-seqread.h"

typedef struct /* information stored for each node of the lcp interval tree */
{
  unsigned long *filenumdist;
} ShulengthdistDfsinfo;

typedef struct  /* global information */
{
  unsigned long numofdbfiles;
  const GtEncseq *encseq;
  unsigned long **shulengthdist;
} Shulengthdiststate;

#include "esa-dfs.h"

static void resetfilenumdist(ShulengthdistDfsinfo *father,
                             unsigned long numofdbfiles)
{
  unsigned long idx;

  for (idx = 0; idx < numofdbfiles; idx++)
  {
    father->filenumdist[idx] = 0;
  }
}

static Dfsinfo *allocateDfsinfo(Dfsstate *astate)
{
  ShulengthdistDfsinfo *dfsinfo;
  Shulengthdiststate *state = (Shulengthdiststate*) astate;

  dfsinfo = gt_malloc(sizeof(*dfsinfo));
  dfsinfo->filenumdist
    = gt_malloc(sizeof(*dfsinfo->filenumdist) * state->numofdbfiles);
  resetfilenumdist(dfsinfo,state->numofdbfiles);
  return (Dfsinfo*) dfsinfo;
}

static void freeDfsinfo(Dfsinfo *adfsinfo, GT_UNUSED Dfsstate *state)
{
  ShulengthdistDfsinfo *dfsinfo = (ShulengthdistDfsinfo*) adfsinfo;;
  gt_free(dfsinfo->filenumdist);
  gt_free(dfsinfo);
}

static int processleafedge(bool firstsucc,
                           unsigned long fatherdepth,
                           Dfsinfo *afather,
                           unsigned long leafnumber,
                           Dfsstate *astate,
                           GT_UNUSED GtError *err)
{
  Shulengthdiststate *state = (Shulengthdiststate*) astate;
  ShulengthdistDfsinfo *father = (ShulengthdistDfsinfo*) afather;
#undef SKDEBUG
#ifdef SKDEBUG
  printf("processleafedge %lu firstsucc=%s, "
         " depth(father)= %lu\n",
         leafnumber,
         firstsucc ? "true" : "false",
         fatherdepth);
#endif
  if (fatherdepth > 0)
  {
    unsigned long filenum;

    filenum = gt_encseq_filenum(state->encseq,leafnumber);
    if (firstsucc)
    {
      resetfilenumdist(father,state->numofdbfiles);
    } else
    {
      unsigned long idx;

      for (idx = 0; idx < state->numofdbfiles; idx++)
      {
        if (idx != filenum)
        {
          state->shulengthdist[idx][filenum] += fatherdepth + 1;
        }
      }
    }
    father->filenumdist[filenum]++;
  } else
  {
    resetfilenumdist(father,state->numofdbfiles);
  }
  return 0;
}

static int processbranchedge(bool firstsucc,
                             unsigned long fatherdepth,
                             Dfsinfo *afather,
                             Dfsinfo *ason,
                             Dfsstate *astate,
                             GT_UNUSED GtError *err)
{
  Shulengthdiststate *state = (Shulengthdiststate*) astate;
  ShulengthdistDfsinfo *father = (ShulengthdistDfsinfo*) afather;

#ifdef SKDEBUG
  printf("processbranchedge firstsucc=%s, depth(father)=%lu\n",
         firstsucc ? "true" : "false",fatherdepth);
#endif
  if (fatherdepth > 0)
  {
    ShulengthdistDfsinfo *son = (ShulengthdistDfsinfo*) ason;
    unsigned long idx;

    if (firstsucc)
    {
      resetfilenumdist(father,state->numofdbfiles);
    } else
    {
      unsigned long idx1, idx2;

      for (idx1=0; idx1 < state->numofdbfiles; idx1++)
      {
        for (idx2=0; idx2 < state->numofdbfiles; idx2++)
        {
          if (idx1 != idx2)
          {
            state->shulengthdist[idx1][idx2]
              += (fatherdepth + 1) * son->filenumdist[idx2];
          }
        }
      }
      for (idx = 0; idx < state->numofdbfiles; idx++)
      {
        father->filenumdist[idx] += son->filenumdist[idx];
      }
    }
  } else
  {
    resetfilenumdist(father,state->numofdbfiles);
  }
  return 0;
}

int gt_sastream2shulengthdist(Sequentialsuffixarrayreader *ssar,
                              const GtEncseq *encseq,
                              GtLogger *logger,
                              GtError *err)
{
  Shulengthdiststate *state;
  bool haserr = false;
  unsigned long idx1, idx2;

  state = gt_malloc(sizeof(*state));
  state->numofdbfiles = gt_encseq_num_of_files(encseq);
  state->encseq = encseq;
  gt_array2dim_malloc(state->shulengthdist,state->numofdbfiles,
                      state->numofdbfiles);
  for (idx1=0; idx1 < state->numofdbfiles; idx1++)
  {
    for (idx2=0; idx2 < state->numofdbfiles; idx2++)
    {
      state->shulengthdist[idx1][idx2] = 0;
    }
  }
  if (gt_depthfirstesa(ssar,
                       allocateDfsinfo,
                       freeDfsinfo,
                       processleafedge,
                       processbranchedge,
                       NULL,
                       NULL,
                       NULL,
                       (Dfsstate*) state,
                       logger,
                       err) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    for (idx1=0; idx1 < state->numofdbfiles; idx1++)
    {
      for (idx2=0; idx2 < state->numofdbfiles; idx2++)
      {
        if (idx1 != idx2)
        {
          printf("%lu %lu %lu\n",idx1,idx2,state->shulengthdist[idx1][idx2]);
        }
      }
    }
  }
  gt_array2dim_delete(state->shulengthdist);
  gt_free(state);
  return haserr ? -1 : 0;
}
