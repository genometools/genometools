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
#include "core/seqiterator_sequence_buffer.h"
#include "esa-seqread.h"
#include "esa-splititv.h"
#define SKDEBUG
#ifdef SKDEBUG
#include "core/encseq.h"
#endif

typedef struct /* information stored for each node of the lcp interval tree */
{
  unsigned long *filenumdist;
  unsigned long id;
} ShulengthdistDfsinfo;

typedef struct  /* global information */
{
  unsigned long numofdbfiles,
                lastleafnumber,
                **shulengthdist;
  const GtEncseq *encseq;
  unsigned long nextid;
} Shulengthdiststate;

#include "esa-dfs.h"

static void resetfilenumdist(ShulengthdistDfsinfo *father,
                             unsigned long numofdbfiles)
{
  unsigned long idx;

  printf("reset node %lu\n",father->id);
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
  dfsinfo->id = state->nextid++;
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

static void contribute(unsigned long **shulengthdist,unsigned long idx1,
                       unsigned long idx2,unsigned long value)
{
#ifdef SKDEBUG
  printf("add[%lu][%lu]+=%lu\n",idx1,idx2,value);
#endif
  shulengthdist[idx1][idx2] += value;
}

static void shownode(const Shulengthdiststate *state,
                     const char *kind,
                     const ShulengthdistDfsinfo *node)
{
  unsigned long idx;

  printf("%s(id=%lu):",kind,node->id);
  for (idx=0; idx < state->numofdbfiles; idx++)
  {
    if (node->filenumdist[idx] > 0)
    {
      printf(" %lu->%lu",idx,node->filenumdist[idx]);
    }
  }
  printf("\n");
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
#ifdef SKDEBUG
  printf("processleafedge %lu firstsucc=%s, "
         " depth(father)= %lu, path=",
         leafnumber,
         firstsucc ? "true" : "false",
         fatherdepth);
  if (fatherdepth > 0)
  {
    gt_encseq_showatstartposwithdepth(stdout,
                                      state->encseq,
                                      GT_READMODE_FORWARD,
                                      leafnumber,
                                      fatherdepth);
  }
  printf("\n");
#endif
  if (fatherdepth > 0)
  {
    unsigned long filenum;

    filenum = gt_encseq_filenum(state->encseq,leafnumber);
    if (firstsucc)
    {
      resetfilenumdist(father,state->numofdbfiles);
      shownode(state,"father",father);
    } else
    {
      unsigned long idx;

      shownode(state,"father",father);
      for (idx = 0; idx < state->numofdbfiles; idx++)
      {
        if (idx != filenum)
        {
          contribute(state->shulengthdist,idx,filenum,fatherdepth + 1);
          if (father->filenumdist[idx] > 0)
          {
            contribute(state->shulengthdist,filenum,idx,
                                            father->filenumdist[idx] *
                                            (fatherdepth + 1));
          }
        }
      }
    }
    father->filenumdist[filenum]++;
  } else
  {
    resetfilenumdist(father,state->numofdbfiles);
  }
  state->lastleafnumber = leafnumber;
  return 0;
}

static void cartproduct(Shulengthdiststate *state,
                        unsigned long depth,
                        ShulengthdistDfsinfo *node1,
                        ShulengthdistDfsinfo *node2)
{
  unsigned long idx1, idx2;

  for (idx1=0; idx1 < state->numofdbfiles; idx1++)
  {
    if (node1->filenumdist[idx1] > 0 && node2->filenumdist[idx1] == 0)
    {
      for (idx2=0; idx2 < state->numofdbfiles; idx2++)
      {
        if (node2->filenumdist[idx2] > 0)
        {
          gt_assert(idx1 != idx2);
          contribute(state->shulengthdist,idx1,idx2,
                     (depth + 1) * node2->filenumdist[idx2]);
        }
      }
    }
  }
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
  printf("processbranchedge firstsucc=%s, depth(father)=%lu,path=",
         firstsucc ? "true" : "false",fatherdepth);
  if (fatherdepth > 0)
  {
    gt_encseq_showatstartposwithdepth(stdout,
                                      state->encseq,
                                      GT_READMODE_FORWARD,
                                      state->lastleafnumber,
                                      fatherdepth);
  }
  printf("\n");
#endif
  if (fatherdepth > 0)
  {
    ShulengthdistDfsinfo *son = (ShulengthdistDfsinfo*) ason;
    unsigned long idx;

    if (firstsucc)
    {
      resetfilenumdist(father,state->numofdbfiles);
      shownode(state,"father",father);
      shownode(state,"son",son);
    } else
    {
      shownode(state,"father",father);
      shownode(state,"son",son);
      cartproduct(state, fatherdepth, father, son);
      cartproduct(state, fatherdepth, son, father);
      for (idx = 0; idx < state->numofdbfiles; idx++)
      {
        father->filenumdist[idx] += son->filenumdist[idx];
      }
      resetfilenumdist(son,state->numofdbfiles);
    }
  } else
  {
    resetfilenumdist(father,state->numofdbfiles);
  }
  return 0;
}

int gt_multiesa2shulengthdist(Sequentialsuffixarrayreader *ssar,
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
  state->nextid = 0;
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

static unsigned long gt_esa2shulengthatposition(const Suffixarray *suffixarray,
                                                unsigned long totallength,
                                                unsigned long offset,
                                                unsigned long left,
                                                unsigned long right,
                                                const GtUchar *qstart,
                                                const GtUchar *qend)
{
  Simplelcpinterval itv;
  const GtUchar *qptr;

  gt_assert(left < right);
  itv.left = left;
  itv.right = right;
  /*printf("\n");*/
  for (qptr = qstart; /* Nothing */; qptr++, offset++)
  {
    if (itv.left <= itv.right)
    {
      /*
      if (qptr < qend)
      {
        printf("read %u\n",(unsigned int) *qptr);
      }
      */
      if (qptr >= qend || ISSPECIAL(*qptr) ||
          !gt_lcpintervalfindcharchildintv(suffixarray->encseq,
                                           suffixarray->readmode,
                                           totallength,
                                           suffixarray->suftab,
                                           &itv,
                                           *qptr,
                                           offset,
                                           itv.left,
                                           itv.right))
      {
        break;
      }
    } else
    {
      break;
    }
  }
  /*printf("add %lu\n",offset+1); */
  return offset+1;
}

static unsigned long gt_esa2shulengthquery(const Suffixarray *suffixarray,
                                           const GtUchar *query,
                                           unsigned long querylen)
{
  const GtUchar *qptr;
  unsigned long totalgmatchlength = 0, gmatchlength, remaining;
  unsigned long totallength = gt_encseq_total_length(suffixarray->encseq);

  for (qptr = query, remaining = querylen; remaining > 0; qptr++, remaining--)
  {
    gmatchlength = gt_esa2shulengthatposition(suffixarray,
                                              totallength,
                                              0,
                                              0,
                                              totallength,
                                              qptr,
                                              query+querylen);
    if (gmatchlength > 0)
    {
      totalgmatchlength += gmatchlength;
    }
  }
  return totalgmatchlength;
}

int gt_esa2shulengthqueryfiles(unsigned long *totalgmatchlength,
                               const Suffixarray *suffixarray,
                               const GtStrArray *queryfilenames,
                               GtError *err)
{
  bool haserr = false;
  GtSeqIterator *seqit;
  const GtUchar *query;
  unsigned long querylen;
  char *desc = NULL;
  int retval;
  GtAlphabet *alphabet;

  gt_error_check(err);
  alphabet = gt_encseq_alphabet(suffixarray->encseq);
  gt_assert(gt_str_array_size(queryfilenames) == 1UL);
  seqit = gt_seqiterator_sequence_buffer_new(queryfilenames, err);
  if (!seqit)
  {
    haserr = true;
  }
  if (!haserr)
  {
    gt_seqiterator_set_symbolmap(seqit, gt_alphabet_symbolmap(alphabet));
    for (; /* Nothing */; )
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
      *totalgmatchlength += gt_esa2shulengthquery(suffixarray,query,querylen);
      gt_free(desc);
    }
    gt_seqiterator_delete(seqit);
  }
  return haserr ? -1 : 0;
}
