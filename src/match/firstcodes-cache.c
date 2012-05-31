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

#include "core/arraydef.h"
#include "core/divmodmul.h"
#include "core/log_api.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "firstcodes-cache.h"
#include "firstcodes-spacelog.h"
#include <limits.h>
#include <math.h>

typedef struct
{
  unsigned long afcindex, code;
} GtIndexwithcode;

struct GtArrayGtIndexwithcode
{
  GtIndexwithcode *spaceGtIndexwithcode;
  unsigned long width, nextfreeGtIndexwithcode, allocatedGtIndexwithcode;
  unsigned int depth;
};

GtArrayGtIndexwithcode *gt_firstcodes_binsearchcache_new(
                                      unsigned long differentcodes,
                                      unsigned int addbscache_depth,
                                      GtFirstcodesspacelog *fcsl)
{
  unsigned int depth = addbscache_depth +
                       (unsigned int) log10((double) differentcodes);
  while (true)
  {
    unsigned long allocated = 1UL << (depth+1);
    if (allocated < differentcodes && differentcodes/allocated >= 32UL)
    {
      GtArrayGtIndexwithcode *binsearchcache;
      size_t allocbytes;

      binsearchcache = gt_malloc(sizeof (*binsearchcache));
      binsearchcache->depth = depth;
      binsearchcache->nextfreeGtIndexwithcode = 0;
      binsearchcache->allocatedGtIndexwithcode = allocated;
      binsearchcache->width = differentcodes/allocated;
      allocbytes = sizeof (*binsearchcache->spaceGtIndexwithcode)
                   * binsearchcache->allocatedGtIndexwithcode;
      gt_assert(binsearchcache->width > 0);
      binsearchcache->spaceGtIndexwithcode = gt_malloc(allocbytes);
      gt_log_log("binsearchcache->depth=%u => %lu bytes",
                 binsearchcache->depth,
                 (unsigned long) allocbytes);
      GT_FCI_ADDWORKSPACE(fcsl,"binsearchcache",allocbytes);
      return binsearchcache;
    } else
    {
      if (depth > 0)
      {
        depth--;
      } else
      {
        break;
      }
    }
  }
  return NULL;
}

unsigned long gt_firstcodes_binsearchcache_width(const GtArrayGtIndexwithcode
                                                 *binsearchcache)
{
  return binsearchcache->width;
}

void gt_firstcodes_binsearchcache_set_index_code(GtArrayGtIndexwithcode
                                                 *binsearchcache,
                                                 unsigned long afcindex,
                                                 unsigned long code)
{
  if (binsearchcache->nextfreeGtIndexwithcode <
      binsearchcache->allocatedGtIndexwithcode)
  {
    binsearchcache->spaceGtIndexwithcode
      [binsearchcache->nextfreeGtIndexwithcode].afcindex = afcindex;
    binsearchcache->spaceGtIndexwithcode
      [binsearchcache->nextfreeGtIndexwithcode++].code = code;
  }
}

void gt_firstcodes_binsearchcache_check(GtArrayGtIndexwithcode *binsearchcache,
                                       GT_UNUSED const unsigned long
                                                                 *allfirstcodes,
                                       unsigned long differentcodes)
{
  if (binsearchcache != NULL)
  {
    unsigned long idx, current = binsearchcache->width;

    gt_assert(binsearchcache->spaceGtIndexwithcode != NULL);
    for (idx=0; idx < binsearchcache->nextfreeGtIndexwithcode &&
                current < differentcodes; idx++)
    {
      gt_assert(binsearchcache->spaceGtIndexwithcode[idx].afcindex == current);
      gt_assert(binsearchcache->spaceGtIndexwithcode[idx].code
                == allfirstcodes[current]);
      current += binsearchcache->width;
    }
  }
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

unsigned long gt_firstcodes_find_accu(unsigned long *foundcode,
                                      const unsigned long *differences,
                                      unsigned long allfirstcodes0,
                                      unsigned long differentcodes,
                                      unsigned long differencemask,
                                      const GtArrayGtIndexwithcode
                                         *binsearchcache,
                                      unsigned long code)
{
  unsigned long leftptr = ULONG_MAX, rightptr = ULONG_MAX,
                foundindex = ULONG_MAX, previouscode = ULONG_MAX;

  if (code <= allfirstcodes0)
  {
    *foundcode = allfirstcodes0;
    return 0;
  }
  *foundcode = ULONG_MAX;
  if (binsearchcache != NULL)
  {
    const GtIndexwithcode *leftic, *midic, *rightic;
    unsigned int depth;

    gt_assert(binsearchcache->spaceGtIndexwithcode != NULL);
    leftic = binsearchcache->spaceGtIndexwithcode;
    rightic = binsearchcache->spaceGtIndexwithcode +
              binsearchcache->nextfreeGtIndexwithcode - 1;
    for (depth = 0; /* Nothing */; depth++)
    {
      midic = leftic + GT_DIV2((unsigned long) (rightic-leftic));
      if (code < midic->code)
      {
        foundindex = midic->afcindex;
        *foundcode = midic->code;
        if (depth < binsearchcache->depth)
        {
          rightic = midic - 1;
        } else
        {
          gt_assert(leftic->afcindex != ULONG_MAX &&
                    rightic->afcindex != ULONG_MAX);
          if (leftic > binsearchcache->spaceGtIndexwithcode)
          {
            leftptr = (leftic-1)->afcindex + 1;
            previouscode = (leftic-1)->code;
          } else
          {
            gt_assert(code > allfirstcodes0);
            leftptr = 1UL;
            previouscode = allfirstcodes0;
          }
          gt_assert(rightic->afcindex > 0);
          rightptr = rightic->afcindex - 1;
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
            gt_assert(leftic->afcindex != ULONG_MAX &&
                      rightic->afcindex != ULONG_MAX);
            leftptr = leftic->afcindex + 1;
            previouscode = leftic->code;
            if (rightic < binsearchcache->spaceGtIndexwithcode +
                          binsearchcache->nextfreeGtIndexwithcode - 1)
            {
              gt_assert((rightic+1)->afcindex > 0);
              rightptr = (rightic+1)->afcindex - 1;
            } else
            {
              rightptr = differentcodes - 1;
            }
            break;
          }
        } else
        {
          gt_assert(midic->afcindex != ULONG_MAX);
          *foundcode = midic->code;
          return midic->afcindex;
        }
      }
    }
    gt_assert(leftptr != ULONG_MAX && rightptr != ULONG_MAX);
  } else
  {
    leftptr = 1UL;
    previouscode = allfirstcodes0;
    rightptr = differentcodes - 1;
  }
  if (leftptr <= rightptr)
  {
    unsigned long idx;

    for (idx = leftptr; idx <= rightptr; idx++)
    {
      previouscode += (differences[idx] & differencemask); /* extract diff */
      if (code <= previouscode)
      {
        *foundcode = previouscode;
        foundindex = idx;
        break;
      }
    }
  }
  return foundindex;
}
