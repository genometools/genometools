/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include <stdio.h>
#include <limits.h>
#include "core/unused_api.h"
#include "core/logger.h"
#include "core/encseq.h"
#include "core/ma_api.h"
#include "sarr-def.h"
#include "emimergeesa.h"
#include "merger-trie.h"
#include "lcpoverflow.h"

#include "esa-map.h"

static void fillandinsert(Mergertrierep *trierep,
                          unsigned int idx,
                          unsigned long suftabvalue,
                          Mergertrienode *node,
                          GT_UNUSED uint64_t ident)
{
  Suffixinfo sinfo;

  sinfo.idx = idx;
  sinfo.startpos = suftabvalue;
#ifdef WITHTRIEIDENT
  sinfo.ident = ident;
#endif
  gt_mergertrie_insertsuffix(trierep,node,&sinfo);
}

static int inputthesequences(unsigned int *numofchars,
                             unsigned long *nextpostable,
                             Suffixarray *suffixarraytable,
                             const GtStrArray *indexnametab,
                             unsigned int demand,
                             GtLogger *logger,
                             GtError *err)
{
  unsigned long idx;
  const char *indexname;

  gt_error_check(err);
  for (idx=0; idx<gt_str_array_size(indexnametab); idx++)
  {
    indexname = gt_str_array_get(indexnametab,idx);
    if (streamsuffixarray(&suffixarraytable[idx],
                          demand,
                          indexname,
                          logger,
                          err) != 0)
    {
      return -1;
    }
    if (idx == 0)
    {
      *numofchars =
            gt_alphabet_num_of_chars(
                     gt_encseq_alphabet(suffixarraytable[idx].encseq));
    }
    nextpostable[idx] = 0;
  }
  return 0;
}

static int insertfirstsuffixes(Mergertrierep *trierep,
                               unsigned long *nextpostable,
                               Suffixarray *suffixarraytable,
                               unsigned int numofindexes,
                               GtError *err)
{
  unsigned int idx;
  unsigned long suftabvalue;
  int retval;

  gt_error_check(err);
  for (idx=0; idx<numofindexes; idx++)
  {
    retval
      = gt_readnextfromstream_GtUlong(&suftabvalue,
                                  &suffixarraytable[idx].suftabstream_GtUlong);
    if (retval == 0)
    {
      gt_error_set(err,"file %s: line %d: unexpected end of file when "
                        "reading suftab",__FILE__,__LINE__);
      return -2;
    }
    nextpostable[idx]++;
    fillandinsert(trierep,
                  idx,
                  suftabvalue,
                  trierep->root,
                  (uint64_t) idx);
  }
  return 0;
}

/*@null@*/ static Mergertrienode *findlargestnodeleqlcpvalue(
                                          Mergertrienode *smallest,
                                          unsigned long lcpvalue,
                                          GtError *err)
{
  Mergertrienode *tmp;

  gt_error_check(err);
  for (tmp = smallest->parent; tmp != NULL; tmp = tmp->parent)
  {
    if (tmp->depth <= lcpvalue)
    {
      return tmp;
    }
  }
  gt_error_set(err,"path does not contain a node of depth <= %lu",
                lcpvalue);
  return NULL;
}

int gt_emissionmergedesa_stepdeleteandinsertothersuffixes(
                                 Emissionmergedesa *emmesa, GtError *err)
{
  Mergertrienode *tmpsmallestleaf, *tmplcpnode;
  Largelcpvalue tmpexception;
  GtUchar tmpsmalllcpvalue;
  int retval;
  unsigned long tmpsuftabvalue,
         tmplcpvalue,
         tmplastbranchdepth;
  unsigned int tmpidx;

  gt_error_check(err);
  for (emmesa->buf.nextstoreidx = 0;
      emmesa->numofentries > 0 &&
      emmesa->buf.nextstoreidx < (unsigned int) SIZEOFMERGERESULTBUFFER;
      emmesa->buf.nextstoreidx++)
  {
    tmpsmallestleaf = gt_mergertrie_findsmallestnode(&emmesa->trierep);
    tmplastbranchdepth = tmpsmallestleaf->parent->depth;
    tmpidx = tmpsmallestleaf->suffixinfo.idx;
    emmesa->buf.suftabstore[emmesa->buf.nextstoreidx].idx = tmpidx;
    emmesa->buf.suftabstore[emmesa->buf.nextstoreidx].startpos
      = tmpsmallestleaf->suffixinfo.startpos;
    if (emmesa->nextpostable[tmpidx] >
       gt_encseq_total_length(emmesa->suffixarraytable[tmpidx].encseq))
    {
      gt_mergertrie_deletesmallestpath(tmpsmallestleaf,&emmesa->trierep);
      emmesa->numofentries--;
    } else
    {
      retval = gt_readnextfromstream_GtUchar(&tmpsmalllcpvalue,
                                       &emmesa->suffixarraytable[tmpidx].
                                                lcptabstream);
      if (retval < 0)
      {
        return -1;
      }
      if (retval == 0)
      {
        gt_error_set(err,"file %s: line %d: unexpected end of file when "
                        "reading lcptab",__FILE__,__LINE__);
        return -2;
      }
      if (tmpsmalllcpvalue == LCPOVERFLOW)
      {
        retval = gt_readnextfromstream_Largelcpvalue(
                               &tmpexception,
                               &emmesa->suffixarraytable[tmpidx].llvtabstream);
        if (retval < 0)
        {
          return -3;
        }
        if (retval == 0)
        {
          gt_error_set(err,"file %s: line %d: unexpected end of file when "
                        "reading llvtab",__FILE__,__LINE__);
          return -4;
        }
        tmplcpvalue = tmpexception.value;
      } else
      {
        tmplcpvalue = (unsigned long) tmpsmalllcpvalue;
      }
      if (tmplcpvalue > tmplastbranchdepth)
      {
        tmplastbranchdepth = tmplcpvalue;
      }
      tmplcpnode = findlargestnodeleqlcpvalue(tmpsmallestleaf,tmplcpvalue,err);
      retval = gt_readnextfromstream_GtUlong(&tmpsuftabvalue,
                                         &emmesa->suffixarraytable[tmpidx].
                                         suftabstream_GtUlong);
      if (retval == 0)
      {
        gt_error_set(err,"file %s: line %d: unexpected end of file when "
                      "reading suftab",__FILE__,__LINE__);
        return -6;
      }
      emmesa->nextpostable[tmpidx]++;
      fillandinsert(&emmesa->trierep,
                    tmpidx,
                    tmpsuftabvalue,
                    tmplcpnode,
                    emmesa->ident++);
      tmpsmallestleaf = gt_mergertrie_findsmallestnode(&emmesa->trierep);
      gt_mergertrie_deletesmallestpath(tmpsmallestleaf,&emmesa->trierep);
    }
    if (emmesa->numofentries > 0)
    {
      emmesa->buf.lcptabstore[emmesa->buf.nextstoreidx]
        = tmplastbranchdepth;
      emmesa->buf.lastpage = false;
    } else
    {
      emmesa->buf.lastpage = true;
    }
  }
  return 0;
}

int gt_emissionmergedesa_init(Emissionmergedesa *emmesa,
                           const GtStrArray *indexnametab,
                           unsigned int demand,
                           GtLogger *logger,
                           GtError *err)
{
  unsigned int numofindexes;
  bool haserr = false;

  numofindexes = (unsigned int) gt_str_array_size(indexnametab);
  emmesa->buf.nextaccessidx = emmesa->buf.nextstoreidx = 0;
  emmesa->numofindexes = numofindexes;
  emmesa->numofentries = numofindexes;
  emmesa->ident = (uint64_t) numofindexes;
  emmesa->trierep.encseqreadinfo = NULL;
  emmesa->suffixarraytable = gt_malloc(sizeof *emmesa->suffixarraytable
                                       * numofindexes);
  emmesa->nextpostable = gt_malloc(sizeof *emmesa->nextpostable * numofindexes);
  if (inputthesequences(&emmesa->numofchars,
                        emmesa->nextpostable,
                        emmesa->suffixarraytable,
                        indexnametab,
                        demand,
                        logger,
                        err) != 0)
  {
    haserr = true;
    return -1;
  }
  if (!haserr && numofindexes > 1U)
  {
    unsigned int idx;

    emmesa->trierep.encseqreadinfo
      = gt_malloc(sizeof *emmesa->trierep.encseqreadinfo * numofindexes);
    for (idx = 0; idx < numofindexes; idx++)
    {
      emmesa->trierep.encseqreadinfo[idx].encseqptr
        = emmesa->suffixarraytable[idx].encseq;
      emmesa->trierep.encseqreadinfo[idx].readmode
        = emmesa->suffixarraytable[idx].readmode;
    }
    gt_mergertrie_initnodetable(&emmesa->trierep,(unsigned long) numofindexes,
                             numofindexes);
    if (insertfirstsuffixes(&emmesa->trierep,
                           emmesa->nextpostable,
                           emmesa->suffixarraytable,
                           numofindexes,
                           err) != 0)
    {
      gt_free(emmesa->trierep.encseqreadinfo);
      emmesa->trierep.encseqreadinfo = NULL;
      haserr = true;
    }
  }
  if (haserr)
  {
    gt_free(emmesa->suffixarraytable);
    gt_free(emmesa->nextpostable);
  }
  return haserr ? -1 : 0;
}

void gt_emissionmergedesa_wrap(Emissionmergedesa *emmesa)
{
  unsigned int idx;

  for (idx = 0; idx < emmesa->numofindexes; idx++)
  {
    gt_freesuffixarray(emmesa->suffixarraytable + idx);
  }
  gt_free(emmesa->suffixarraytable);
  gt_free(emmesa->trierep.encseqreadinfo);
  emmesa->trierep.encseqreadinfo = NULL;
  if (emmesa->numofindexes > 1U)
  {
    gt_mergertrie_delete(&emmesa->trierep);
  }
  gt_free(emmesa->nextpostable);
}
