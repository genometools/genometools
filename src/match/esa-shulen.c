/*
  Copyright (c) 2010 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
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
#include "core/format64.h"
#undef SHUDEBUG
#ifdef SHUDEBUG
#include "core/encseq.h"
#endif
#include "esa-seqread.h"
#include "esa-splititv.h"
#include "shu_unitfile.h"
#include "esa-shulen.h"

typedef struct /* information stored for each node of the lcp interval tree */
{
  unsigned long *gnumdist;
#ifdef SHUDEBUG
  unsigned long id;
#endif
} BUinfo_shulen;

typedef struct  /* global information */
{
  unsigned long numofdbfiles,
                lastleafnumber;
  uint64_t **shulengthdist;
  const GtEncseq *encseq;
  unsigned long *file_to_genome_map;
#ifdef SHUDEBUG
  unsigned long nextid;
#endif
} BUstate_shulen;

static void resetgnumdist_shulen(BUinfo_shulen *father,
                                    unsigned long numofdbfiles)
{
  unsigned long idx;

#ifdef SHUDEBUG
  printf("reset node %lu\n",father->id);
#endif
  for (idx = 0; idx < numofdbfiles; idx++)
  {
    father->gnumdist[idx] = 0;
  }
}

static void initBUinfo_shulen(BUinfo_shulen *buinfo,
                              GT_UNUSED BUstate_shulen *state)
{
#ifdef SHUDEBUG
  buinfo->id = state->nextid++;
#endif
  buinfo->gnumdist = NULL;
}

static void freeBUinfo_shulen(BUinfo_shulen *buinfo,
                              GT_UNUSED BUstate_shulen *state)
{
  gt_free(buinfo->gnumdist);
}

static void contribute_shulen(GT_UNUSED int line,
                              uint64_t **shulengthdist,
                              unsigned long referidx,
                              unsigned long shulenidx,
                              unsigned long count,
                              unsigned long depth)
{
#ifdef SHUDEBUG
  printf("line %d: add[%lu][%lu]+=count=%lu*depth=%lu\n",line,referidx,
                                                        shulenidx,count,depth);
#endif
  shulengthdist[referidx][shulenidx] += count * depth;
}

#ifdef SHUDEBUG
static void shownode(int line,const BUstate_shulen *state,
                     const char *kind,
                     const BUinfo_shulen *node)
{
  unsigned long idx;

  printf("line %d: %s(id=%lu,numofdbfiles=%lu):",line,kind,node->id,
                                                 state->numofdbfiles);
  for (idx=0; idx < state->numofdbfiles; idx++)
  {
    if (node->gnumdist[idx] > 0)
    {
      printf(" %lu->%lu",idx,node->gnumdist[idx]);
    }
  }
  printf("\n");
}
#endif

static int processleafedge_shulen(bool firstsucc,
                                  unsigned long fatherdepth,
                                  GT_UNUSED unsigned long fatherlb,
                                  BUinfo_shulen *father,
                                  unsigned long leafnumber,
                                  BUstate_shulen *state,
                                  GT_UNUSED GtError *err)
{
  unsigned long idx, gnum;

#ifdef SHUDEBUG
  printf("processleafedge %lu firstsucc=%s, "
         " depth(father)=%lu, path=",
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
  if (state->file_to_genome_map != NULL)
  {
    gnum = state->file_to_genome_map[gt_encseq_filenum(state->encseq,
                                                       leafnumber)];
  } else
  {
    gnum = gt_encseq_filenum(state->encseq,leafnumber);
  }
  if (firstsucc)
  {
    gt_assert(father != NULL);
    if (father->gnumdist == NULL)
    {
      father->gnumdist
        = gt_malloc(sizeof (*father->gnumdist) * state->numofdbfiles);
    }
    resetgnumdist_shulen(father,state->numofdbfiles);
#ifdef SHUDEBUG
    shownode(__LINE__,state,"father",father);
#endif
  } else
  {
#ifdef SHUDEBUG
    shownode(__LINE__,state,"father",father);
#endif
    for (idx = 0; idx < state->numofdbfiles; idx++)
    {
      if (idx != gnum)
      {
        if (father->gnumdist[idx] > 0)
        {
          contribute_shulen(__LINE__,state->shulengthdist,idx,gnum,
                            1UL,fatherdepth+1);
          if (father->gnumdist[gnum] == 0)
          {
            contribute_shulen(__LINE__,state->shulengthdist,gnum,idx,
                              father->gnumdist[idx],fatherdepth + 1);
          }
        }
      }
    }
  }
  father->gnumdist[gnum]++;
#ifdef SHUDEBUG
  printf("gnumdist[id=%lu,filenum=%lu]=%lu\n",father->id,gnum,
                                              father->gnumdist[gnum]);
#endif
  state->lastleafnumber = leafnumber;
  return 0;
}

static void cartproduct_shulen(BUstate_shulen *state,
                               unsigned long depth,
                               const BUinfo_shulen *node1,
                               const BUinfo_shulen *node2)
{
  unsigned long referidx, shulenidx;

  for (referidx=0; referidx < state->numofdbfiles; referidx++)
  {
    if (node1->gnumdist[referidx] > 0 && node2->gnumdist[referidx] == 0)
    {
      for (shulenidx=0; shulenidx < state->numofdbfiles; shulenidx++)
      {
        if (node2->gnumdist[shulenidx] > 0)
        {
          gt_assert(referidx != shulenidx);
          contribute_shulen(__LINE__,state->shulengthdist,referidx,shulenidx,
                            node2->gnumdist[shulenidx],depth + 1);
        }
      }
    }
  }
}

static int processbranchingedge_shulen(bool firstsucc,
                                       unsigned long fatherdepth,
                                       GT_UNUSED unsigned long fatherlp,
                                       BUinfo_shulen *father,
                                       GT_UNUSED unsigned long sondepth,
                                       GT_UNUSED unsigned long sonlb,
                                       GT_UNUSED unsigned long sonrb,
                                       BUinfo_shulen *son,
                                       BUstate_shulen *state,
                                       GT_UNUSED GtError *err)
{
  unsigned long idx;

#ifdef SHUDEBUG
  printf("%s firstsucc=%s, depth(father)=%lu,path=",__func__,
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
  if (firstsucc)
  {
    gt_assert(father != NULL);
    if (father->gnumdist == NULL)
    {
      father->gnumdist
        = gt_malloc(sizeof (*father->gnumdist) * state->numofdbfiles);
      resetgnumdist_shulen(father,state->numofdbfiles);
    }
#ifdef SHUDEBUG
    shownode(__LINE__,state,"father",father);
#endif
  } else
  {
#ifdef SHUDEBUG
    gt_assert(father != NULL);
    shownode(__LINE__,state,"father",father);
    gt_assert(son != NULL);
    shownode(__LINE__,state,"son",son);
#endif
    cartproduct_shulen(state, fatherdepth, father, son);
    cartproduct_shulen(state, fatherdepth, son, father);
  }
  if (son != NULL)
  {
    for (idx = 0; idx < state->numofdbfiles; idx++)
    {
      father->gnumdist[idx] += son->gnumdist[idx];
      son->gnumdist[idx] = 0;
#ifdef SHUDEBUG
      printf("gnumdist[id=%lu,filenum=%lu]=%lu\n",father->id,idx,
                                                father->gnumdist[idx]);
      printf("gnumdist[id=%lu,filenum=%lu]=0\n",son->id,idx);
#endif
    }
  }
  return 0;
}

static int processcompletenode_shulen(GT_UNUSED unsigned long lcp,
                                      GT_UNUSED unsigned long lb,
                                      GT_UNUSED unsigned long rb,
                                      GT_UNUSED BUinfo_shulen *info,
                                      GT_UNUSED BUstate_shulen *state,
                                      GT_UNUSED GtError *err)
{
  return 0;
}

#include "esa-bottomup-shulen.inc"

int gt_multiesa2shulengthdist(Sequentialsuffixarrayreader *ssar,
                              const GtEncseq *encseq,
                              GtError *err)
{
  BUstate_shulen *state;
  bool haserr = false;
  unsigned long referidx, shulenidx;

  state = gt_malloc(sizeof (*state));
  state->numofdbfiles = gt_encseq_num_of_files(encseq);
  state->encseq = encseq;
#ifdef SHUDEBUG
  state->nextid = 0;
#endif
  gt_array2dim_malloc(state->shulengthdist,state->numofdbfiles,
                      state->numofdbfiles);
  for (referidx=0; referidx < state->numofdbfiles; referidx++)
  {
    for (shulenidx=0; shulenidx < state->numofdbfiles; shulenidx++)
    {
      state->shulengthdist[referidx][shulenidx] = 0;
    }
  }
  if (gt_esa_bottomup_shulen(ssar, state, err) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    printf("%lu\n",state->numofdbfiles);
    for (shulenidx=0; shulenidx < state->numofdbfiles; shulenidx++)
    {
      printf("%lu\t",shulenidx);
      for (referidx=0; referidx < state->numofdbfiles; referidx++)
      {
        if (referidx != shulenidx)
        {
          printf(Formatuint64_t"\t",
                 PRINTuint64_tcast(state->shulengthdist[referidx][shulenidx]));
        } else
        {
          printf("0\t");
        }
      }
      printf("\n");
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
    if (ISSPECIAL(*qptr))
    {
      gmatchlength = 0;
    } else
    {
      gmatchlength = gt_esa2shulengthatposition(suffixarray,
                                              totallength,
                                              0,
                                              0,
                                              totallength,
                                              qptr,
                                              query+querylen);
    }
    totalgmatchlength += gmatchlength;
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
    }
    gt_seqiterator_delete(seqit);
  }
  return haserr ? -1 : 0;
}

int gt_get_multiesashulengthdist(Sequentialsuffixarrayreader *ssar,
                                 const GtEncseq *encseq,
                                 uint64_t **shulen,
                                 struct GtShuUnitFileInfo_tag *unit_info,
                                 GtError *err)
{
  BUstate_shulen *state;
  bool haserr = false;

  state = gt_malloc(sizeof (*state));
  state->numofdbfiles = unit_info->num_of_genomes;
  state->file_to_genome_map = unit_info->map_files;
  state->encseq = encseq;
#ifdef SHUDEBUG
  state->nextid = 0;
#endif
  state->shulengthdist = shulen;
  if (gt_esa_bottomup_shulen(ssar, state, err) != 0)
  {
    haserr = true;
  }
  gt_free(state);
  return haserr ? -1 : 0;
}
