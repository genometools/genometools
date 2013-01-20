/*
  Copyright (c) 2009-2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2009-2011 Center for Bioinformatics, University of Hamburg

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

#include "core/arraydef.h"
#include "core/assert_api.h"
#include "match/esa-seqread.h"
#include "core/encseq.h"
#include "core/progressbar.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "core/intbits.h"
#include "match/esa-splititv.h"
#include "match/rdj-revcompl-def.h"
#include "match/rdj-ovlfind-gusfield.h"

typedef struct
{
  unsigned long offset, leftbound,
                rightbound, specials_leftbound;
  bool visited;
} RdjGusfieldIndexBounds;

GT_DECLAREARRAYSTRUCT(RdjGusfieldIndexBounds);

static inline
void maximize_child_offset(RdjGusfieldIndexBounds *child,
                           const Suffixarray *sarr)
{
  unsigned long idx, lcp;
  child->offset = lcptable_get(sarr, child->leftbound+1);
  for (idx = child->leftbound+2; idx <= child->rightbound; idx++)
  {
    lcp = lcptable_get(sarr, idx);
    if (lcp < child->offset)
    {
      child->offset = lcp;
    }
  }
}

static inline
void rdj_gusfield_processleafedge(unsigned long leafnumber,
                                  const GtEncseq *encseq,
                                  unsigned long nofsequences,
                                  GtArrayGtUlong *matchlen_stacks,
                                  GtBitsequence *has_match,
                                  GtBitsequence *sspbittab,
                                  bool find_submaximal,
                                  unsigned long firstrevcompl,
                                  GtSpmproc proc,
                                  void* procdata)
{
  unsigned long prefix_seqnum, suffix_seqnum, corrected_suffix_seqnum;
  bool suffixseq_direct = true, prefixseq_direct = true;
  unsigned long *matchlen;
  unsigned long stack_pos;

  if (leafnumber == 0 || (bool)GT_ISIBITSET(sspbittab, leafnumber-1))
  {
    prefix_seqnum = gt_encseq_seqnum(encseq, leafnumber);
    if (firstrevcompl != 0)
    {
      GT_READJOINER_CORRECT_REVCOMPL(prefix_seqnum, firstrevcompl,
        nofsequences, prefixseq_direct);
    }
    for (suffix_seqnum = 0; suffix_seqnum < nofsequences; suffix_seqnum++)
    {
      if ((bool)GT_ISIBITSET(has_match, suffix_seqnum))
      {
        corrected_suffix_seqnum = suffix_seqnum;
        if (firstrevcompl != 0)
        {
          GT_READJOINER_CORRECT_REVCOMPL(corrected_suffix_seqnum,
            firstrevcompl, nofsequences, suffixseq_direct);
          if (!GT_READJOINER_IS_CORRECT_REVCOMPL_CASE(corrected_suffix_seqnum,
            suffixseq_direct, prefix_seqnum, prefixseq_direct))
            continue;
        }
        if (!find_submaximal)
        {
          matchlen = matchlen_stacks[suffix_seqnum].spaceGtUlong +
                     matchlen_stacks[suffix_seqnum].nextfreeGtUlong - 1;
          if (suffix_seqnum == prefix_seqnum)
          {
            if (*matchlen == gt_encseq_seqlength(encseq, suffix_seqnum))
            {
              /* read overlaps with whole itself */
              if (matchlen_stacks[suffix_seqnum].nextfreeGtUlong == 1UL)
              {
                /* no other overlap */
                continue;
              }
              else
              {
                matchlen = matchlen_stacks[suffix_seqnum].spaceGtUlong +
                           matchlen_stacks[suffix_seqnum].nextfreeGtUlong - 2;
              }
            }
          }
          proc(corrected_suffix_seqnum, prefix_seqnum,
               (unsigned long) *matchlen, suffixseq_direct,
               prefixseq_direct, procdata);
        }
        else /* process whole stack, not only the top */
        {
          for (stack_pos = 0;
               stack_pos < matchlen_stacks[suffix_seqnum].nextfreeGtUlong;
               stack_pos++)
          {
            matchlen = matchlen_stacks[suffix_seqnum].spaceGtUlong + stack_pos;
            if (suffix_seqnum == prefix_seqnum)
            {
              if (*matchlen == gt_encseq_seqlength(encseq, suffix_seqnum))
                /* skip ovl of read with whole itself */
                continue;
            }
            proc(corrected_suffix_seqnum, prefix_seqnum,
                 (unsigned long) *matchlen, suffixseq_direct,
                 prefixseq_direct, procdata);
          }
        }
      }
    }
  }
}

/* search terminal edges and push/pop offset value on/from spm stacks */
static inline
void processspecials(unsigned long specials_leftbound,
                     unsigned long rightbound,
                     unsigned long offset,
                     bool visited,
                     const unsigned long *suftab,
                     const GtEncseq *encseq,
                     unsigned long totallength,
                     unsigned long nofsequences,
                     GtArrayGtUlong *matchlen_stacks,
                     GtBitsequence *has_match,
                     GtBitsequence *sspbittab)
{
  unsigned long pos, leafnumber;
  unsigned long seqnum;

  gt_assert(offset>0);

  for (pos = specials_leftbound; pos <= rightbound; pos++)
  {
    leafnumber = suftab[pos];
    if ((bool)GT_ISIBITSET(sspbittab, leafnumber + offset))
    {
      seqnum = (leafnumber == totallength)
               ? nofsequences-1
               : gt_encseq_seqnum(encseq,leafnumber);
      if (!visited)
      {
        GT_SETIBIT(has_match, seqnum);
        GT_STOREINARRAY(&(matchlen_stacks[seqnum]),
                        GtUlong, 1, offset);
      }
      else
      {
        gt_assert(matchlen_stacks[seqnum].nextfreeGtUlong>0);
        matchlen_stacks[seqnum].nextfreeGtUlong--;
        if (matchlen_stacks[seqnum].nextfreeGtUlong == 0)
        {
          GT_UNSETIBIT(has_match, seqnum);
        }
      }
    }
  }
}

void gt_rdj_gusfield(Sequentialsuffixarrayreader *ssar,
    unsigned long min_length, bool find_submaximal, bool show_progressbar,
    unsigned long firstrevcompl, GtSpmproc proc, void* procdata)
{

  /* index */
  const Suffixarray *sarr;
  const GtEncseq *encseq;
  const unsigned long *suftab;
  unsigned long nofsequences;
  GtReadmode readmode;
  unsigned long totallength;
  unsigned int numofchars;
  GtBitsequence *sspbittab;

  /* splitting */
  GtArrayRdjGusfieldIndexBounds stack;
  GtArrayBoundswithchar bwci;
  RdjGusfieldIndexBounds parent, child, *parent_stack_ptr;
  unsigned long seqpos_idx;
  unsigned long idx;

  /* stacks for spm determination */
  GtArrayGtUlong *matchlen_stacks;
  GtBitsequence *has_match;

  /* progressbar */
  unsigned long long progress;

  gt_assert(min_length>0);
  gt_assert(ssar != NULL);

  sarr             = gt_suffixarraySequentialsuffixarrayreader(ssar);
  encseq           = gt_encseqSequentialsuffixarrayreader(ssar);
  suftab           = gt_suftabSequentialsuffixarrayreader(ssar);
  nofsequences     = gt_encseq_num_of_sequences(encseq);
  readmode         = gt_readmodeSequentialsuffixarrayreader(ssar);
  totallength      = gt_encseq_total_length(encseq);
  numofchars       = gt_alphabet_num_of_chars(gt_encseq_alphabet(encseq));

  /* prepare sspbittab */
  GT_INITBITTAB(sspbittab, totallength+1);
  for (idx = 1UL; idx <= nofsequences - 1; idx++)
  {
    GT_SETIBIT(sspbittab, gt_encseq_seqstartpos(encseq, idx) - 1);
  }
  GT_SETIBIT(sspbittab, totallength);

  matchlen_stacks = gt_malloc(sizeof *matchlen_stacks * nofsequences);
  GT_INITBITTAB(has_match, nofsequences);

  for (idx = 0UL; idx < nofsequences; idx++)
  {
    GT_INITARRAY(&(matchlen_stacks[idx]), GtUlong);
  }

  if (show_progressbar)
  {
    progress = 0ULL;
    gt_progressbar_start(&progress, (unsigned long long)totallength);
  }

  GT_INITARRAY(&stack, RdjGusfieldIndexBounds);

  /* prepare bwci array */
  bwci.spaceBoundswithchar = gt_malloc(sizeof *bwci.spaceBoundswithchar *
                                       (numofchars + 1));
  bwci.nextfreeBoundswithchar = 0UL;
  bwci.allocatedBoundswithchar = (unsigned long) (numofchars + 1);

  /* get rid of compiler warnings */
  parent.specials_leftbound = 0UL;
  child.specials_leftbound = 0UL;

  /* start at [0,n] */
  parent.offset = 0UL;
  parent.leftbound = 0UL;
  parent.rightbound = totallength;
  parent.visited = false;
  GT_STOREINARRAY(&stack, RdjGusfieldIndexBounds, 1UL, parent);

  while (stack.nextfreeRdjGusfieldIndexBounds > 0UL)
  {
    parent_stack_ptr = stack.spaceRdjGusfieldIndexBounds +
                          stack.nextfreeRdjGusfieldIndexBounds - 1;
    parent = *parent_stack_ptr;

    if (parent.visited)
    {
      /* search terminal edges and pop spm stacks */
      if (parent.specials_leftbound <= parent.rightbound &&
           parent.offset >= min_length)
        processspecials(parent.specials_leftbound, parent.rightbound,
                        parent.offset, parent.visited, suftab, encseq,
                        totallength, nofsequences, matchlen_stacks, has_match,
                        sspbittab);
      stack.nextfreeRdjGusfieldIndexBounds--;
      continue;
    }

    parent_stack_ptr->visited = true;

    /* split interval */
    gt_lcpintervalsplitwithoutspecial(&bwci, encseq, readmode,
                                   totallength, suftab,
                                   parent.offset,
                                   parent.leftbound,
                                   parent.rightbound);

    if (bwci.nextfreeBoundswithchar > 0)
    {
      parent.specials_leftbound = bwci.spaceBoundswithchar[
                                    bwci.nextfreeBoundswithchar-1].rbound+1;
    }
    else
    {
      parent.specials_leftbound = parent.leftbound;
    }
    parent_stack_ptr->specials_leftbound = parent.specials_leftbound;

    /* search terminal edges and push on spm stacks */
    if (parent.specials_leftbound <= parent.rightbound &&
         parent.offset >= min_length)
      processspecials(parent.specials_leftbound, parent.rightbound,
                      parent.offset, parent.visited, suftab, encseq,
                      totallength, nofsequences, matchlen_stacks, has_match,
                      sspbittab);

    for (idx = bwci.nextfreeBoundswithchar; idx > 0; idx--)
    {
      child.leftbound = bwci.spaceBoundswithchar[idx-1].lbound;
      child.rightbound = bwci.spaceBoundswithchar[idx-1].rbound;
      child.visited = false;

      if (child.leftbound == child.rightbound)
      {
        if (parent.offset >= min_length)
          rdj_gusfield_processleafedge(suftab[child.leftbound], encseq,
                          nofsequences, matchlen_stacks, has_match, sspbittab,
                          find_submaximal, firstrevcompl, proc, procdata);
        if (show_progressbar) progress++;
      }
      else
      {
        maximize_child_offset(&child, sarr);
        GT_STOREINARRAY(&stack, RdjGusfieldIndexBounds, 1, child);
      }
    }

    for (seqpos_idx = parent.specials_leftbound;
         seqpos_idx <= parent.rightbound;
         seqpos_idx++)
    {
      if (parent.offset >= min_length)
        rdj_gusfield_processleafedge(suftab[seqpos_idx], encseq,
                        nofsequences, matchlen_stacks, has_match, sspbittab,
                        find_submaximal, firstrevcompl, proc, procdata);
      if (show_progressbar) progress++;
    }
  }

  GT_FREEARRAY(&stack, RdjGusfieldIndexBounds);
  GT_FREEARRAY(&bwci,Boundswithchar);

  if (show_progressbar) gt_progressbar_stop();

  for (idx = 0; idx < nofsequences; idx++)
  {
    GT_FREEARRAY(&(matchlen_stacks[idx]), GtUlong);
  }
  gt_free(has_match);
  gt_free(sspbittab);
  gt_free(matchlen_stacks);
}
