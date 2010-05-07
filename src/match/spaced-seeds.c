/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include <stdarg.h>
#include "core/assert_api.h"
#include "core/types_api.h"
#include "core/unused_api.h"
#include "core/chardef.h"
#include "core/ma_api.h"
#include "core/intbits.h"

#include "absdfstrans-imp.h"

typedef struct
{
  bool pathmatches;
} SpseLimdfsstate;

typedef struct
{
  GtBitsequence seedbitvector;
  unsigned long seedweight;
  const GtUchar *pattern;
} SpseLimdfsconstinfo;

#ifdef SKDEBUG

static void spse_showLimdfsstate(const DECLAREPTRDFSSTATE(aliascol),
                                unsigned long currentdepth,
                                GT_UNUSED const Limdfsconstinfo *mti)
{
  const SpseLimdfsstate *col = (const SpseLimdfsstate *) aliascol;

  printf("at depth %lu (pathmatches=%s)\n",currentdepth,
                                           col->pathmatches ? "true" : "false");
}

#endif

static Limdfsconstinfo *spse_allocatedfsconstinfo(
                               GT_UNUSED unsigned int alphasize)
{
  SpseLimdfsconstinfo *mti = gt_malloc(sizeof (SpseLimdfsconstinfo));
  return (Limdfsconstinfo*) mti;
}

static void spse_initdfsconstinfo(Limdfsconstinfo *mt,
                                  unsigned int alphasize,
                                  ...)
                                 /* Variable argument list is as follows:
                                    const GtUchar *pattern,
                                    Bitsequence seedbitvector,
                                    unsigned long seedweight
                                 */
{
  va_list ap;
  SpseLimdfsconstinfo *mti = (SpseLimdfsconstinfo*) mt;

  va_start(ap,alphasize);
  mti->pattern = va_arg(ap, const GtUchar *);
  mti->seedbitvector = va_arg(ap, GtBitsequence);
  mti->seedweight = va_arg(ap, unsigned long);
  va_end(ap);
}

static void spse_freedfsconstinfo(Limdfsconstinfo **mtiptr)
{
  gt_free(*mtiptr);
  *mtiptr = NULL;
}

static void spse_initLimdfsstate(DECLAREPTRDFSSTATE(aliascolumn),
                                 GT_UNUSED Limdfsconstinfo *mti)
{
  SpseLimdfsstate *column = (SpseLimdfsstate *) aliascolumn;

  column->pathmatches = true;
}

static void spse_fullmatchLimdfsstate(Limdfsresult *limdfsresult,
                                      DECLAREPTRDFSSTATE(aliascolumn),
                                      GT_UNUSED unsigned long leftbound,
                                      GT_UNUSED unsigned long rightbound,
                                      GT_UNUSED unsigned long width,
                                      unsigned long currentdepth,
                                      Limdfsconstinfo *mt)
{
  SpseLimdfsstate *limdfsstate = (SpseLimdfsstate *) aliascolumn;
  SpseLimdfsconstinfo *mti = (SpseLimdfsconstinfo*) mt;

  if (limdfsstate->pathmatches)
  {
    if (currentdepth == mti->seedweight)
    {
      limdfsresult->status = Limdfssuccess;
      limdfsresult->pprefixlen = mti->seedweight;
      limdfsresult->distance = 0;
      return;
    }
    if (currentdepth < mti->seedweight)
    {
      limdfsresult->status = Limdfscontinue;
      return;
    }
  }
  limdfsresult->status = Limdfsstop;
}

static bool setpathmatch(GtBitsequence seedbitvector,
                         const GtUchar *pattern,
                         unsigned long currentdepth,
                         GtUchar currentchar)
{
  return (!GT_ISBITSET(seedbitvector,currentdepth-1) ||
          currentchar == pattern[currentdepth-1]) ? true : false;
}

static void spse_nextLimdfsstate(const Limdfsconstinfo *mt,
                                 DECLAREPTRDFSSTATE(aliasoutcol),
                                 unsigned long currentdepth,
                                 GtUchar currentchar,
                                 GT_UNUSED const DECLAREPTRDFSSTATE(aliasincol))
{
  SpseLimdfsstate *outcol = (SpseLimdfsstate *) aliasoutcol;
#ifndef NDEBUG
  const SpseLimdfsstate *incol = (const SpseLimdfsstate *) aliasincol;
#endif
  SpseLimdfsconstinfo *mti = (SpseLimdfsconstinfo*) mt;

  gt_assert(ISNOTSPECIAL(currentchar));
  gt_assert(currentdepth > 0);
  gt_assert(incol->pathmatches);

  outcol->pathmatches = setpathmatch(mti->seedbitvector,
                                     mti->pattern,
                                     currentdepth,
                                     currentchar);
}

static void spse_inplacenextLimdfsstate(const Limdfsconstinfo *mt,
                                        DECLAREPTRDFSSTATE(aliascol),
                                        unsigned long currentdepth,
                                        GtUchar currentchar)
{
  SpseLimdfsstate *col = (SpseLimdfsstate *) aliascol;
  const SpseLimdfsconstinfo *mti = (const SpseLimdfsconstinfo*) mt;

  gt_assert(ISNOTSPECIAL(currentchar));
  gt_assert(currentdepth > 0);
  col->pathmatches = setpathmatch(mti->seedbitvector,
                                  mti->pattern,
                                  currentdepth,
                                  currentchar);
}

const AbstractDfstransformer *gt_spse_AbstractDfstransformer(void)
{
  static const AbstractDfstransformer spse_adfst =
  {
    sizeof (SpseLimdfsstate),
    spse_allocatedfsconstinfo,
    spse_initdfsconstinfo,
    NULL, /* no extractdfsconstinfo */
    spse_freedfsconstinfo,
    spse_initLimdfsstate,
    NULL,
    NULL,
    NULL,
    spse_fullmatchLimdfsstate,
    spse_nextLimdfsstate,
    spse_inplacenextLimdfsstate,
#ifdef SKDEBUG
    spse_showLimdfsstate,
#endif
  };
  return &spse_adfst;
}
