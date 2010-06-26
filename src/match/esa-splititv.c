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

#include "core/chardef.h"
#include "core/divmodmul.h"
#include "core/types_api.h"

#include "esa-splititv.h"

#define SEQUENCE(ENCSEQ,POS) (((POS) == totallength) \
                             ? (GtUchar) SEPARATOR\
                             : gt_encseq_get_encoded_char(ENCSEQ, \
                                                                 POS, \
                                                                 readmode))

static unsigned long lcpintervalfindrightbound(const GtEncseq *encseq,
                                        GtReadmode readmode,
                                        unsigned long totallength,
                                        const ESASuffixptr *suftab,
                                        GtUchar cc,
                                        unsigned long offset,
                                        unsigned long left,
                                        unsigned long right)
{
  unsigned long pos, mid;
  GtUchar midcc;

  while (right > left+1)
  {
    mid = GT_DIV2(left+right);
    pos = ESASUFFIXPTRGET(suftab,mid) + offset;
    midcc = SEQUENCE(encseq,pos);
    if (cc < midcc)
    {
      right = mid;
    } else
    {
      left = mid;
    }
  }
  return left;
}

bool gt_lcpintervalfindcharchildintv(const GtEncseq *encseq,
                                     GtReadmode readmode,
                                     unsigned long totallength,
                                     const ESASuffixptr *suftab,
                                     Simplelcpinterval *itv,
                                     GtUchar cc,
                                     unsigned long offset,
                                     unsigned long left,
                                     unsigned long right)
{
  GtUchar leftcc, rightcc;
  unsigned long pos, rightbound, leftbound = left;

  pos = ESASUFFIXPTRGET(suftab,right) + offset;
  rightcc = SEQUENCE(encseq,pos);
  while (true)
  {
    pos = ESASUFFIXPTRGET(suftab,leftbound) + offset;
    leftcc = SEQUENCE(encseq,pos);
    if (leftcc == rightcc)
    {
      break;
    }
    rightbound = lcpintervalfindrightbound(encseq,readmode,
                                           totallength,suftab,leftcc,
                                           offset,leftbound,right);
    if (leftcc == cc)
    {
      itv->left = leftbound;
      itv->right = rightbound;
      return true;
    }
    if (leftcc > cc)
    {
      return false;
    }
    leftbound = rightbound+1;
  }
  if (leftcc == cc)
  {
    itv->left = leftbound;
    itv->right = right;
    return true;
  }
  return false;
}

#define ADDCURRENTLBOUND(V)\
        bwci->spaceBoundswithchar[bwci->nextfreeBoundswithchar].lbound = V

#define ADDPREVIOUSRBOUND(V)\
        if (bwci->nextfreeBoundswithchar > 0)\
        {\
          bwci->spaceBoundswithchar[bwci->nextfreeBoundswithchar-1].rbound = V;\
        }

#define ADDCURRENTINCHAR(V)\
        bwci->spaceBoundswithchar[bwci->nextfreeBoundswithchar++].inchar = V

void gt_lcpintervalsplitwithoutspecial(GtArrayBoundswithchar *bwci,
                                       const GtEncseq *encseq,
                                       GtReadmode readmode,
                                       unsigned long totallength,
                                       const ESASuffixptr *suftab,
                                       unsigned long parentoffset,
                                       unsigned long parentleft,
                                       unsigned long parentright)
{
  GtUchar leftcc, rightcc;
  unsigned long rightbound = 0, leftbound = parentleft;

  /* call gt_lcpintervalextendlcp and verify if interval can be extended by
     some character */
  bwci->nextfreeBoundswithchar = 0;
  rightcc = SEQUENCE(encseq,ESASUFFIXPTRGET(suftab,parentright) + parentoffset);
  while (true)
  {
    leftcc = SEQUENCE(encseq,ESASUFFIXPTRGET(suftab,leftbound) + parentoffset);
    gt_assert(bwci->nextfreeBoundswithchar < bwci->allocatedBoundswithchar);
    if (ISSPECIAL(leftcc))
    {
      ADDPREVIOUSRBOUND(rightbound);
      ADDCURRENTLBOUND(rightbound+1);
      return;
    }
    ADDPREVIOUSRBOUND(leftbound-1);
    ADDCURRENTLBOUND(leftbound);
    ADDCURRENTINCHAR(leftcc);
    if (leftcc == rightcc)
    {
      break;
    }
    rightbound = lcpintervalfindrightbound(encseq,readmode,totallength,suftab,
                                           leftcc,parentoffset,
                                           leftbound,parentright);
    leftbound = rightbound+1;
  }
  gt_assert(bwci->nextfreeBoundswithchar < bwci->allocatedBoundswithchar);
  ADDPREVIOUSRBOUND(parentright);
  ADDCURRENTLBOUND(parentright+1);
}

GtUchar gt_lcpintervalextendlcp(const GtEncseq *encseq,
                                GtReadmode readmode,
                                const ESASuffixptr *suftab,
                                unsigned long totallength,
                                GtUchar alphasize,
                                unsigned long parentoffset,
                                unsigned long parentleft,
                                unsigned long parentright)
{
  GtUchar ccl, ccr;

  ccl = SEQUENCE(encseq,ESASUFFIXPTRGET(suftab,parentleft) + parentoffset);
  ccr = SEQUENCE(encseq,ESASUFFIXPTRGET(suftab,parentright) + parentoffset);
  if (ccl != ccr || ISSPECIAL(ccl))
  {
    return alphasize;
  }
  gt_assert(ccl < alphasize);
  return ccl;
}
