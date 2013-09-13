/*
  Copyright (c) 2013 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

static GtUword lcpintervalfindfirst(const GtEncseq *encseq,
                                          GtReadmode readmode,
                                          GtUword totallength,
                                          const ESASuffixptr *suftab,
                                          GtUchar cc,
                                          GtUword offset,
                                          GtUword left,
                                          GtUword right)
{
  GtUword found = ULONG_MAX;

  while (left <= right)
  {
    GtUword mid = left + GT_DIV2(right - left + 1);
    GtUword pos = ESASUFFIXPTRGET(suftab,mid) + offset;
    GtUchar midcc = SEQUENCE(encseq,pos);
    if (midcc < cc)
    {
      left = mid + 1;
    } else
    {
      if (cc == midcc)
      {
        found = mid;
      }
      if (mid == 0)
      {
        break;
      }
      right = mid - 1;
    }
  }
  return found;
}

static GtUword lcpintervalfindlast(const GtEncseq *encseq,
                                         GtReadmode readmode,
                                         GtUword totallength,
                                         const ESASuffixptr *suftab,
                                         GtUchar cc,
                                         GtUword offset,
                                         GtUword left,
                                         GtUword right)
{
  GtUword found = ULONG_MAX;

  while (left <= right)
  {
    GtUword mid = left + GT_DIV2(right - left + 1);
    GtUword pos = ESASUFFIXPTRGET(suftab,mid) + offset;
    GtUchar midcc = SEQUENCE(encseq,pos);
    if (cc < midcc)
    {
      if (mid == 0)
      {
        break;
      }
      right = mid - 1;
    } else
    {
      if (cc == midcc)
      {
        found = mid;
      }
      left = mid + 1;
    }
  }
  return found;
}

bool gt_lcpintervalfindcharchildintv_simple(const GtEncseq *encseq,
                                     GtReadmode readmode,
                                     GtUword totallength,
                                     const ESASuffixptr *suftab,
                                     Simplelcpinterval *itv,
                                     GtUchar cc,
                                     GtUword offset,
                                     GtUword left,
                                     GtUword right)
{
  GtUword pos;

  pos = ESASUFFIXPTRGET(suftab,left) + offset;
  if (cc < SEQUENCE(encseq,pos))
  {
    return false;
  }
  pos = ESASUFFIXPTRGET(suftab,right) + offset;
  if (cc > SEQUENCE(encseq,pos))
  {
    return false;
  }
  itv->left = lcpintervalfindfirst(encseq,
                                   readmode,
                                   totallength,
                                   suftab,
                                   cc,
                                   offset,
                                   left,
                                   right);
  if (itv->left == ULONG_MAX)
  {
    return false;
  }
  itv->right = lcpintervalfindlast(encseq,
                                   readmode,
                                   totallength,
                                   suftab,
                                   cc,
                                   offset,
                                   itv->left + 1,
                                   right);
  if (itv->right == ULONG_MAX)
  {
    itv->right = itv->left;
  }
  return true;
}

bool gt_lcpintervalfindcharchildintv_withcheck(const GtEncseq *encseq,
                                     GtReadmode readmode,
                                     GtUword totallength,
                                     const ESASuffixptr *suftab,
                                     Simplelcpinterval *itv,
                                     GtUchar cc,
                                     GtUword offset,
                                     GtUword left,
                                     GtUword right)
{
  Simplelcpinterval itv2;
  bool occurs = gt_lcpintervalfindcharchildintv(encseq,
                                          readmode,
                                      totallength,
                                      suftab,
                                      itv,
                                      cc,
                                      offset,left,right);
  bool occurs2 = gt_lcpintervalfindcharchildintv_simple(encseq,
                                                  readmode,
                                              totallength,
                                              suftab,
                                              &itv2,
                                              cc,
                                              offset,left,right);
  if ((occurs && !occurs2) || (!occurs && occurs2))
  {
    fprintf(stderr,"occurs = %s, occurs2 = %s\n",occurs ? "true" : "false",
                                                 occurs2 ? "true" : "false");
    exit(EXIT_FAILURE);
  }
  if (occurs && itv->left != itv2.left)
  {
    fprintf(stderr,"left = " GT_WU " != " GT_WU " = left2\n",
            itv->left,itv2.left);
    exit(EXIT_FAILURE);
  }
  if (occurs && itv->right != itv2.right)
  {
    fprintf(stderr,"right = " GT_WU " != " GT_WU " = right2\n",
            itv->right,itv2.right);
    exit(EXIT_FAILURE);
  }
  return occurs;
}
