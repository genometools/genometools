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
#include "core/seq_iterator_sequence_buffer_api.h"
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
} GtBUinfo_shulen;

struct GtBUstate_shulen /* global information */
{
  unsigned long numofdbfiles;
  uint64_t **shulengthdist;
  const GtEncseq *encseq;
  unsigned long *file_to_genome_map;
#undef GENOMEDIFF_PAPER_IMPL
#ifdef GENOMEDIFF_PAPER_IMPL
  unsigned long *leafdist;
#endif
#ifdef SHUDEBUG
  unsigned long lastleafnumber,
                nextid;
#endif
  /* the remaining components are used for the genomediff function
     called from the suffixerator */
  unsigned long idxoffset,
                previousbucketlastsuffix;
  bool firstedgefromroot;
  GtShuUnitFileInfo *unit_info;
  void *stack;
};

static void resetgnumdist_shulen(GtBUinfo_shulen *father,
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

static void initBUinfo_shulen(GtBUinfo_shulen *buinfo,
                              GT_UNUSED GtBUstate_shulen *state)
{
#ifdef SHUDEBUG
  buinfo->id = state->nextid++;
#endif
  buinfo->gnumdist = NULL;
}

static void freeBUinfo_shulen(GtBUinfo_shulen *buinfo,
                              GT_UNUSED GtBUstate_shulen *state)
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
static void shownode(int line,
                     const GtBUstate_shulen *state,
                     const char *kind,
                     const GtBUinfo_shulen *node)
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

static void cartproduct_shulen(GtBUstate_shulen *state,
                               unsigned long depth,
                               const unsigned long *refnumdist,
                               const unsigned long *querynumdist)
{
  unsigned long referidx, shulenidx;

  for (referidx=0; referidx < state->numofdbfiles; referidx++)
  {
    if (refnumdist[referidx] > 0 && querynumdist[referidx] == 0)
    {
      for (shulenidx=0; shulenidx < state->numofdbfiles; shulenidx++)
      {
        if (querynumdist[shulenidx] > 0)
        {
          gt_assert(referidx != shulenidx);
          contribute_shulen(__LINE__,
                            state->shulengthdist,
                            referidx,
                            shulenidx,
                            querynumdist[shulenidx],
                            depth + 1);
        }
      }
    }
  }
}

static void shu_compute_leaf_edge_contrib(GtBUstate_shulen *state,
                                          const unsigned long *fathernumdist,
                                          unsigned long gnum,
                                          unsigned long fatherdepth)
{
  unsigned long idx;
#ifdef GENOMEDIFF_PAPER_IMPL
  gt_assert(state->leafdist != NULL);
  for (idx = 0; idx < state->numofdbfiles; idx++)
  {
    state->leafdist[idx] = 0;
  }
  state->leafdist[gnum] = 1UL;
  cartproduct_shulen(state,fatherdepth,fathernumdist,state->leafdist);
  cartproduct_shulen(state,fatherdepth,state->leafdist,fathernumdist);
#else
  for (idx = 0; idx < state->numofdbfiles; idx++)
  {
    if (idx != gnum && fathernumdist[idx] > 0)
    {
      contribute_shulen(__LINE__,
                        state->shulengthdist,
                        idx,
                        gnum,
                        1UL,
                        fatherdepth+1);
      if (fathernumdist[gnum] == 0)
      {
        contribute_shulen(__LINE__,
                          state->shulengthdist,
                          gnum,
                          idx,
                          fathernumdist[idx],
                          fatherdepth + 1);
      }
    }
  }
#endif
}

static int processleafedge_shulen(bool firstsucc,
                                  unsigned long fatherdepth,
                                  GtBUinfo_shulen *father,
                                  unsigned long leafnumber,
                                  GtBUstate_shulen *state,
                                  GT_UNUSED GtError *err)
{
  unsigned long gnum;

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
    shu_compute_leaf_edge_contrib(state,father->gnumdist,gnum,fatherdepth);
  }
  father->gnumdist[gnum]++;
#ifdef SHUDEBUG
  printf("gnumdist[id=%lu,filenum=%lu]=%lu\n",father->id,gnum,
                                              father->gnumdist[gnum]);
  state->lastleafnumber = leafnumber;
#endif
  return 0;
}

static int processbranchingedge_shulen(bool firstsucc,
                                       unsigned long fatherdepth,
                                       GtBUinfo_shulen *father,
                                       GT_UNUSED unsigned long sondepth,
                                       GT_UNUSED unsigned long sonwidth,
                                       GtBUinfo_shulen *son,
                                       GtBUstate_shulen *state,
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
    cartproduct_shulen(state, fatherdepth, father->gnumdist, son->gnumdist);
    cartproduct_shulen(state, fatherdepth, son->gnumdist, father->gnumdist);
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

static uint64_t **shulengthdist_new(unsigned long numofdbfiles)
{
  unsigned long idx1, idx2;
  uint64_t **shulengthdist;

  gt_array2dim_malloc(shulengthdist,numofdbfiles,numofdbfiles);
  for (idx1=0; idx1 < numofdbfiles; idx1++)
  {
    for (idx2=0; idx2 < numofdbfiles; idx2++)
    {
      shulengthdist[idx1][idx2] = 0;
    }
  }
  return shulengthdist;
}

/*
  XXX: check what kind of const annotation is necesseay to
  make the contents of shulengthdist non-writable.
*/

static void shulengthdist_print(const GtStrArray *file_names,
                                const uint64_t * const*shulengthdist,
                                unsigned long numofdbfiles)
{
  unsigned long idx1, idx2;

  /*shulengthdist[0][0] = 0; for testing */
  printf("# sum of shulen\n%lu\n",numofdbfiles);
  for (idx2=0; idx2 < numofdbfiles; idx2++)
  {
    if (file_names != NULL)
    {
      printf("%s\t",gt_str_array_get(file_names,idx2));
    } else
    {
      printf("%lu\t",idx2);
    }
    for (idx1=0; idx1 < numofdbfiles; idx1++)
    {
      if (idx1 != idx2)
      {
        printf(Formatuint64_t"\t",
               PRINTuint64_tcast(shulengthdist[idx1][idx2]));
      } else
      {
        printf("0.000000\t");
      }
    }
    printf("\n");
  }
}

#include "esa-bottomup-shulen.inc"

int gt_multiesa2shulengthdist_print(Sequentialsuffixarrayreader *ssar,
                                    const GtEncseq *encseq,
                                    GtError *err)
{
  GtBUstate_shulen *state;
  bool haserr = false;

  state = gt_malloc(sizeof (*state));
  state->numofdbfiles = gt_encseq_num_of_files(encseq);
  state->encseq = encseq;
#ifdef GENOMEDIFF_PAPER_IMPL
  state->leafdist = gt_malloc(sizeof (*state->leafdist) * state->numofdbfiles);
#endif
#ifdef SHUDEBUG
  state->nextid = 0;
#endif
  state->shulengthdist = shulengthdist_new(state->numofdbfiles);
  if (gt_esa_bottomup_shulen(ssar, state, err) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    shulengthdist_print(NULL,(const uint64_t * const*) state->shulengthdist,
                        state->numofdbfiles);
  }
  gt_array2dim_delete(state->shulengthdist);
#ifdef GENOMEDIFF_PAPER_IMPL
  gt_free(state->leafdist);
#endif
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
  seqit = gt_seq_iterator_sequence_buffer_new(queryfilenames, err);
  if (!seqit)
  {
    haserr = true;
  }
  if (!haserr)
  {
    gt_seq_iterator_set_symbolmap(seqit, gt_alphabet_symbolmap(alphabet));
    for (; /* Nothing */; )
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
      *totalgmatchlength += gt_esa2shulengthquery(suffixarray,query,querylen);
    }
    gt_seq_iterator_delete(seqit);
  }
  return haserr ? -1 : 0;
}

int gt_multiesa2shulengthdist(Sequentialsuffixarrayreader *ssar,
                              const GtEncseq *encseq,
                              uint64_t **shulen,
                              const GtShuUnitFileInfo *unit_info,
                              GtError *err)
{
  GtBUstate_shulen *bustate;
  bool haserr = false;

  bustate = gt_malloc(sizeof (*bustate));
  bustate->numofdbfiles = unit_info->num_of_genomes;
  bustate->file_to_genome_map = unit_info->map_files;
  bustate->encseq = encseq;
#ifdef GENOMEDIFF_PAPER_IMPL
  bustate->leafdist
    = gt_malloc(sizeof (*bustate->leafdist) * bustate->numofdbfiles);
#endif
#ifdef SHUDEBUG
  bustate->nextid = 0;
#endif
  bustate->shulengthdist = shulen;
  if (gt_esa_bottomup_shulen(ssar, bustate, err) != 0)
  {
    haserr = true;
  }
#ifdef GENOMEDIFF_PAPER_IMPL
  gt_free(bustate->leafdist);
#endif
  gt_free(bustate);
  return haserr ? -1 : 0;
}

GtBUstate_shulen *gt_sfx_multiesashulengthdist_new(const GtEncseq *encseq,
                                            GenomediffInfo *gd_info)
{
  GtBUstate_shulen *bustate;

  bustate = gt_malloc(sizeof (*bustate));
  bustate->encseq = encseq;
  bustate->previousbucketlastsuffix = ULONG_MAX;
  bustate->idxoffset = 0;
  bustate->firstedgefromroot = false;
#ifdef SHUDEBUG
  bustate->nextid = 0;
#endif
  if (gd_info == NULL)
    bustate->unit_info = gt_shu_unit_info_new(encseq);
  else
    bustate->unit_info = gd_info->unit_info;

  bustate->numofdbfiles = gt_encseq_num_of_files(encseq);
#ifdef GENOMEDIFF_PAPER_IMPL
  bustate->leafdist
    = gt_malloc(sizeof (*bustate->leafdist) * bustate->numofdbfiles);
#endif
  bustate->file_to_genome_map = bustate->unit_info->map_files;
  if (gd_info == NULL)
    bustate->shulengthdist = shulengthdist_new(bustate->numofdbfiles);
  else
    bustate->shulengthdist = gd_info->shulensums;

  bustate->stack = (void *) gt_GtArrayGtBUItvinfo_new_shulen();
  return bustate;
}

#include "esa-bottomup-shulen-RAM.inc"

int gt_sfx_multiesa2shulengthdist(GtBUstate_shulen *bustate,
                                  const unsigned long *bucketofsuffixes,
                                  const uint32_t *bucketofsuffixes_uint32,
                                  const GtLcpvaluetype *lcptab_bucket,
                                  unsigned long numberofsuffixes,
                                  GtError *err)
{
  bool haserr = false;

  if (bustate->previousbucketlastsuffix != ULONG_MAX &&
      gt_esa_bottomup_RAM_previousfromlast_shulen(
                                 bustate->previousbucketlastsuffix,
                                 (unsigned long) lcptab_bucket[0],
                                 bustate->stack,
                                 bustate,
                                 err) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (gt_esa_bottomup_RAM_shulen(bucketofsuffixes,
                                   bucketofsuffixes_uint32,
                                   lcptab_bucket,
                                   numberofsuffixes,
                                   bustate->stack,
                                   bustate,
                                   err) != 0)
    {
      haserr = true;
    }
  }
  bustate->idxoffset += numberofsuffixes;
  return haserr ? -1 : 0;
}

int gt_sfx_multiesa2shulengthdist_last(GtBUstate_shulen *bustate,GtError *err)
{
  if (bustate->previousbucketlastsuffix != ULONG_MAX &&
      gt_esa_bottomup_RAM_previousfromlast_shulen(
                                 bustate->previousbucketlastsuffix,
                                 0,
                                 bustate->stack,
                                 bustate,
                                 err) != 0)
  {
    return -1;
  }
  return 0;
}

void gt_sfx_multiesashulengthdist_delete(GtBUstate_shulen *bustate,
                                         GenomediffInfo *gd_info)
{
  if (bustate == NULL)
  {
    return;
  }
  gt_assert(bustate->shulengthdist != NULL);
  if (gd_info == NULL) {
    shulengthdist_print(bustate->unit_info->file_names,
                        (const uint64_t * const*) bustate->shulengthdist,
                        bustate->numofdbfiles);
    gt_array2dim_delete(bustate->shulengthdist);
    gt_shu_unit_info_delete(bustate->unit_info);
  }
  gt_GtArrayGtBUItvinfo_delete_shulen(bustate->stack,bustate);
#ifdef GENOMEDIFF_PAPER_IMPL
  gt_free(bustate->leafdist);
#endif
  gt_free(bustate);
}
