/*
  Copyright (c) 2011-2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011-2012 Center for Bioinformatics, University of Hamburg

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

#include <math.h>
#include <limits.h>
#include "core/ma.h"
#include "core/log_api.h"
#include "core/arraydef.h"
#include "core/divmodmul.h"
#include "firstcodes-spacelog.h"
#include "firstcodes-cache.h"

typedef struct
{
  const unsigned long *afcptr;
  unsigned long code;
} GtIndexwithcode;

struct GtArrayGtIndexwithcode
{
  GtIndexwithcode *spaceGtIndexwithcode;
  unsigned long width, nextfreeGtIndexwithcode, allocatedGtIndexwithcode,
                allfirstcodes0; /* entry at index 0 */
  unsigned int depth;
};

GtArrayGtIndexwithcode *gt_firstcodes_binsearchcache_new(
                                      const unsigned long *allfirstcodes,
                                      unsigned long differentcodes,
                                      unsigned int addbscache_depth,
                                      GtFirstcodesspacelog *fcsl)
{
  size_t allocbytes = 0;
  GtArrayGtIndexwithcode *binsearchcache;

  binsearchcache = gt_malloc(sizeof (*binsearchcache));
  binsearchcache->depth
    = addbscache_depth + (unsigned int) log10((double) differentcodes);
  binsearchcache->nextfreeGtIndexwithcode = 0;
  binsearchcache->allocatedGtIndexwithcode = 1UL << (binsearchcache->depth+1);
  binsearchcache->width
    = differentcodes/binsearchcache->allocatedGtIndexwithcode;
  binsearchcache->allfirstcodes0 = allfirstcodes[0];
  if (binsearchcache->allocatedGtIndexwithcode < differentcodes)
  {
    unsigned long idx, current = binsearchcache->width;

    allocbytes = sizeof (*binsearchcache->spaceGtIndexwithcode)
                 * binsearchcache->allocatedGtIndexwithcode;
    binsearchcache->spaceGtIndexwithcode = gt_malloc(allocbytes);
    for (idx=0; idx < binsearchcache->allocatedGtIndexwithcode; idx++)
    {
      gt_assert(current < differentcodes);
      binsearchcache->spaceGtIndexwithcode
               [binsearchcache->nextfreeGtIndexwithcode].afcptr
                  = allfirstcodes + current;
      binsearchcache->spaceGtIndexwithcode
               [binsearchcache->nextfreeGtIndexwithcode++].code
                  = allfirstcodes[current];
      current += binsearchcache->width;
    }
  } else
  {
    binsearchcache->spaceGtIndexwithcode = NULL;
  }
  gt_log_log("binsearchcache->depth=%u => %lu bytes",
             binsearchcache->depth,
             (unsigned long) allocbytes);
  GT_FCI_ADDWORKSPACE(fcsl,"binsearchcache",allocbytes);
  return binsearchcache;
}

void gt_firstcodes_binsearchcache_delete(GtArrayGtIndexwithcode *binsearchcache,
                                         GtFirstcodesspacelog *fcsl)
{
  if (binsearchcache != NULL)
  {
    if (binsearchcache->spaceGtIndexwithcode != NULL)
    {
      GT_FCI_SUBTRACTWORKSPACE(fcsl,"binsearchcache");
      GT_FREEARRAY(binsearchcache,GtIndexwithcode);
    }
    gt_free (binsearchcache);
  }
}

unsigned long gt_firstcodes_binsearchcache_allfirstcodes0(
                    const GtArrayGtIndexwithcode *binsearchcache)
{
  return binsearchcache->allfirstcodes0;
}

const unsigned long *gt_firstcodes_find_accu(
                                 unsigned long *foundcode,
                                 const unsigned long *allfirstcodes,
                                 unsigned long differentcodes,
                                 const GtArrayGtIndexwithcode *binsearchcache,
                                 unsigned long code)
{
  const unsigned long *found = NULL, *leftptr = NULL, *rightptr = NULL;
  unsigned long previouscode = ULONG_MAX;

  if (code <= binsearchcache->allfirstcodes0)
  {
    *foundcode = binsearchcache->allfirstcodes0;
    return allfirstcodes;
  }
  *foundcode = ULONG_MAX;
  if (binsearchcache->spaceGtIndexwithcode != NULL)
  {
    const GtIndexwithcode *leftic, *midic, *rightic;
    unsigned int depth;

    leftic = binsearchcache->spaceGtIndexwithcode;
    rightic = binsearchcache->spaceGtIndexwithcode +
              binsearchcache->nextfreeGtIndexwithcode - 1;
    for (depth = 0; /* Nothing */; depth++)
    {
      midic = leftic + GT_DIV2((unsigned long) (rightic-leftic));
      if (code < midic->code)
      {
        found = midic->afcptr;
        *foundcode = midic->code;
        if (depth < binsearchcache->depth)
        {
          rightic = midic - 1;
        } else
        {
          gt_assert(leftic->afcptr != NULL && rightic->afcptr != NULL);
          if (leftic > binsearchcache->spaceGtIndexwithcode)
          {
            leftptr = (leftic-1)->afcptr + 1;
            previouscode = (leftic-1)->code;
          } else
          {
            gt_assert(code > binsearchcache->allfirstcodes0);
            leftptr = allfirstcodes + 1;
            previouscode = binsearchcache->allfirstcodes0;
          }
          gt_assert(rightic->afcptr > allfirstcodes);
          rightptr = rightic->afcptr - 1;
          break;
        }
      } else
      {
        if (code > midic->code)
        {
          if (depth < binsearchcache->depth)
          {
            leftic = midic + 1;
          } else
          {
            gt_assert(leftic->afcptr != NULL && rightic->afcptr != NULL);
            leftptr = leftic->afcptr + 1;
            previouscode = leftic->code;
            if (rightic < binsearchcache->spaceGtIndexwithcode +
                          binsearchcache->nextfreeGtIndexwithcode - 1)
            {
              gt_assert((rightic+1)->afcptr > allfirstcodes);
              rightptr = (rightic+1)->afcptr - 1;
            } else
            {
              rightptr = allfirstcodes + differentcodes - 1;
            }
            break;
          }
        } else
        {
          gt_assert(midic->afcptr != NULL);
          *foundcode = midic->code;
          return midic->afcptr;
        }
      }
    }
    gt_assert(leftptr != NULL && rightptr != NULL);
  } else
  {
    leftptr = allfirstcodes + 1;
    previouscode = binsearchcache->allfirstcodes0;
    rightptr = allfirstcodes + differentcodes - 1;
  }
  if (leftptr <= rightptr)
  {
    const unsigned long *diff_ptr;

    for (diff_ptr = leftptr; diff_ptr <= rightptr; diff_ptr++)
    {
      previouscode += *diff_ptr; /* extract diff */
      if (code <= previouscode)
      {
        *foundcode = previouscode;
        found = diff_ptr;
        break;
      }
    }
  }
  return found;
}
