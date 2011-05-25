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
#include "esa-bottomup.h"
#include "shu_unitfile.h"
#include "esa-shulen.h"
#include "stamp.h"

typedef struct /* information stored for each node of the lcp interval tree */
{
  unsigned long *filenumdist;
#ifdef SHUDEBUG
  unsigned long id;
#endif
} ShulengthdistBUinfo;

typedef struct  /* global information */
{
  unsigned long numofdbfiles,
                lastleafnumber,
                firstleafcount2omit,
                currentleafcount;
  uint64_t **shulengthdist;
  const GtEncseq *encseq;
  unsigned long *file_to_genome_map;
#ifdef SHUDEBUG
  unsigned long nextid;
#endif
} Shulengthdiststate;

static void shulen_resetfilenumdist(ShulengthdistBUinfo *father,
                                    unsigned long numofdbfiles)
{
  unsigned long idx;

#ifdef SHUDEBUG
  printf("reset node %lu\n",father->id);
#endif
  for (idx = 0; idx < numofdbfiles; idx++)
  {
    father->filenumdist[idx] = 0;
  }
}

static GtBUinfo *shulen_allocateBUinfo(GT_UNUSED GtBUstate *astate)
{
  ShulengthdistBUinfo *dfsinfo;
#ifdef SHUDEBUG
  Shulengthdiststate *state = (Shulengthdiststate*) astate;
#endif

  dfsinfo = gt_malloc(sizeof (*dfsinfo));
#ifdef SHUDEBUG
  dfsinfo->id = state->nextid++;
#endif
  dfsinfo->filenumdist = NULL;
  return (GtBUinfo*) dfsinfo;
}

static void shulen_freeBUinfo(GtBUinfo *adfsinfo, GT_UNUSED GtBUstate *state)
{
  ShulengthdistBUinfo *dfsinfo = (ShulengthdistBUinfo*) adfsinfo;
  gt_free(dfsinfo->filenumdist);
  gt_free(dfsinfo);
}

static void shulen_contribute(GT_UNUSED int line,
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
static void shownode(int line,const Shulengthdiststate *state,
                     const char *kind,
                     const ShulengthdistBUinfo *node)
{
  unsigned long idx;

  printf("line %d: %s(id=%lu,numofdbfiles=%lu):",line,kind,node->id,
                                                 state->numofdbfiles);
  for (idx=0; idx < state->numofdbfiles; idx++)
  {
    if (node->filenumdist[idx] > 0)
    {
      printf(" %lu->%lu",idx,node->filenumdist[idx]);
    }
  }
  printf("\n");
}
#endif

static int shulen_processleafedge(bool firstsucc,
                                  unsigned long fatherdepth,
                                  GT_UNUSED unsigned long fatherlb,
                                  GtBUinfo *afather,
                                  unsigned long leafnumber,
                                  GtBUstate *astate,
                                  GT_UNUSED GtError *err)
{
  Shulengthdiststate *state = (Shulengthdiststate*) astate;
  ShulengthdistBUinfo *father = (ShulengthdistBUinfo*) afather;
  unsigned long idx, filenum;

  if (state->currentleafcount >= state->firstleafcount2omit)
  {
    return 0;
  }
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
  filenum = gt_encseq_filenum(state->encseq,leafnumber);
  if (state->file_to_genome_map != NULL)
  {
    filenum = state->file_to_genome_map[filenum];
  }
  if (firstsucc)
  {
    gt_assert(father != NULL);
    if (father->filenumdist == NULL)
    {
      father->filenumdist
        = gt_malloc(sizeof (*father->filenumdist) * state->numofdbfiles);
    }
    shulen_resetfilenumdist(father,state->numofdbfiles);
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
      if (idx != filenum)
      {
        if (father->filenumdist[idx] > 0)
        {
          shulen_contribute(__LINE__,state->shulengthdist,idx,filenum,
                            1UL,fatherdepth+1);
          if (father->filenumdist[filenum] == 0)
          {
            shulen_contribute(__LINE__,state->shulengthdist,filenum,idx,
                              father->filenumdist[idx],
                              fatherdepth + 1);
          }
        }
      }
    }
  }
  father->filenumdist[filenum]++;
#ifdef SHUDEBUG
  printf("filenumdist[id=%lu,filenum=%lu]=%lu\n",father->id,filenum,
                                                 father->filenumdist[filenum]);
#endif
  state->currentleafcount++;
  state->lastleafnumber = leafnumber;
  return 0;
}

static void shulen_cartproduct(Shulengthdiststate *state,
                               unsigned long depth,
                               const ShulengthdistBUinfo *node1,
                               const ShulengthdistBUinfo *node2)
{
  unsigned long referidx, shulenidx;

  for (referidx=0; referidx < state->numofdbfiles; referidx++)
  {
    if (node1->filenumdist[referidx] > 0 && node2->filenumdist[referidx] == 0)
    {
      for (shulenidx=0; shulenidx < state->numofdbfiles; shulenidx++)
      {
        if (node2->filenumdist[shulenidx] > 0)
        {
          gt_assert(referidx != shulenidx);
          shulen_contribute(__LINE__,state->shulengthdist,referidx,shulenidx,
                            node2->filenumdist[shulenidx],depth + 1);
        }
      }
    }
  }
}

static int shulen_processbranchedge(bool firstsucc,
                                    unsigned long fatherdepth,
                                    GT_UNUSED unsigned long fatherlp,
                                    GtBUinfo *afather,
                                    GT_UNUSED unsigned long sondepth,
                                    GT_UNUSED unsigned long sonlp,
                                    GtBUinfo *ason,
                                    GtBUstate *astate,
                                    GT_UNUSED GtError *err)
{
  Shulengthdiststate *state = (Shulengthdiststate*) astate;
  ShulengthdistBUinfo *father = (ShulengthdistBUinfo*) afather;
  ShulengthdistBUinfo *son = (ShulengthdistBUinfo*) ason;
  unsigned long idx;

#ifdef SHUDEBUG
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
  if (firstsucc)
  {
    gt_assert(father != NULL);
    if (father->filenumdist == NULL)
    {
      father->filenumdist
        = gt_malloc(sizeof (*father->filenumdist) * state->numofdbfiles);
      shulen_resetfilenumdist(father,state->numofdbfiles);
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
    shulen_cartproduct(state, fatherdepth, father, son);
    shulen_cartproduct(state, fatherdepth, son, father);
  }
  if (son != NULL)
  {
    for (idx = 0; idx < state->numofdbfiles; idx++)
    {
      father->filenumdist[idx] += son->filenumdist[idx];
      son->filenumdist[idx] = 0;
#ifdef SHUDEBUG
      printf("filenumdist[id=%lu,filenum=%lu]=%lu\n",father->id,idx,
                                                father->filenumdist[idx]);
      printf("filenumdist[id=%lu,filenum=%lu]=0\n",son->id,idx);
#endif
    }
  }
  return 0;
}

int gt_multiesa2shulengthdist(Sequentialsuffixarrayreader *ssar,
                              const GtEncseq *encseq,
                              GtError *err)
{
  Shulengthdiststate *state;
  bool haserr = false;
  unsigned long referidx, shulenidx;

  state = gt_malloc(sizeof (*state));
  state->numofdbfiles = gt_encseq_num_of_files(encseq);
  state->encseq = encseq;
#ifdef SHUDEBUG
  state->nextid = 0;
#endif
  state->firstleafcount2omit = gt_encseq_total_length(encseq) -
                               gt_encseq_specialcharacters(encseq);
  state->currentleafcount = 0;
  gt_array2dim_malloc(state->shulengthdist,state->numofdbfiles,
                      state->numofdbfiles);
  for (referidx=0; referidx < state->numofdbfiles; referidx++)
  {
    for (shulenidx=0; shulenidx < state->numofdbfiles; shulenidx++)
    {
      state->shulengthdist[referidx][shulenidx] = 0;
    }
  }
  if (gt_esa_bottomup(ssar,
                      shulen_allocateBUinfo,
                      shulen_freeBUinfo,
                      shulen_processleafedge,
                      shulen_processbranchedge,
                      (GtBUstate*) state,
                      err) != 0)
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
  Shulengthdiststate *state;
  bool haserr = false;

  state = gt_malloc(sizeof (*state));
  state->numofdbfiles = unit_info->num_of_genomes;
  state->file_to_genome_map = unit_info->map_files;
  state->encseq = encseq;
#ifdef SHUDEBUG
  state->nextid = 0;
#endif
  state->firstleafcount2omit = gt_encseq_total_length(encseq) -
                               gt_encseq_specialcharacters(encseq);
  state->currentleafcount = 0;
  state->shulengthdist = shulen;
  if (gt_esa_bottomup(ssar,
                      shulen_allocateBUinfo,
                      shulen_freeBUinfo,
                      shulen_processleafedge,
                      shulen_processbranchedge,
                      (GtBUstate*) state,
                      err) != 0)
  {
    haserr = true;
  }
  gt_free(state);
  return haserr ? -1 : 0;
}
