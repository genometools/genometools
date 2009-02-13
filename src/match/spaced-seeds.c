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
#include "core/symboldef.h"
#include "core/unused_api.h"
#include "core/chardef.h"
#include "core/ma.h"
#include "intbits.h"
#include "seqpos-def.h"
#include "absdfstrans-imp.h"

typedef struct
{
  bool pathmatches;
} Limdfsstate;

typedef struct
{
  Bitstring seedbitvector;
  unsigned long seedweight;
  const Uchar *pattern;
} Matchtaskinfo;

#ifdef SKDEBUG

static void spse_showLimdfsstate(const DECLAREPTRDFSSTATE(aliascol),
                                unsigned long currentdepth,
                                GT_UNUSED const void *dfsconstinfo)
{
  const Limdfsstate *col = (const Limdfsstate *) aliascol;

  printf("at depth %lu (pathmatches=%s)\n",currentdepth,
                                            col->pathmatches
                                              ? "true" : "false");
}

#endif

static void *spse_allocatedfsconstinfo(GT_UNUSED unsigned int alphasize)
{
  Matchtaskinfo *mti = gt_malloc(sizeof(Matchtaskinfo));
  return mti;
}

static void spse_initdfsconstinfo(void *dfsconstinfo,
                                  unsigned int alphasize,
                                  ...)
                                 /* Variable argument list is as follows:
                                    const Uchar *pattern,
                                    Bitstring seedbitvector,
                                    unsigned long seedweight
                                 */
{
  va_list ap;
  Matchtaskinfo *mti = (Matchtaskinfo *) dfsconstinfo;

  va_start(ap,alphasize);
  mti->pattern = va_arg(ap, const Uchar *);
  mti->seedbitvector = va_arg(ap, Bitstring);
  mti->seedweight = va_arg(ap, unsigned long);
  va_end(ap);
}

static void spse_freedfsconstinfo(void **dfsconstinfo)
{
  Matchtaskinfo *mti = (Matchtaskinfo *) *dfsconstinfo;

  gt_free(mti);
  *dfsconstinfo = NULL;
}

static void spse_initLimdfsstate(DECLAREPTRDFSSTATE(aliascolumn),
                                 GT_UNUSED void *dfsconstinfo)
{
  Limdfsstate *column = (Limdfsstate *) aliascolumn;

  column->pathmatches = true;
}

static void spse_fullmatchLimdfsstate(Limdfsresult *limdfsresult,
                                      DECLAREPTRDFSSTATE(aliascolumn),
                                      GT_UNUSED Seqpos leftbound,
                                      GT_UNUSED Seqpos rightbound,
                                      GT_UNUSED Seqpos width,
                                      unsigned long currentdepth,
                                      void *dfsconstinfo)
{
  Limdfsstate *limdfsstate = (Limdfsstate *) aliascolumn;
  Matchtaskinfo *mti = (Matchtaskinfo *) dfsconstinfo;

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

static bool setpathmatch(Bitstring seedbitvector,
                         const Uchar *pattern,
                         unsigned long currentdepth,
                         Uchar currentchar)
{
  return (!ISBITSET(seedbitvector,currentdepth-1) ||
          currentchar == pattern[currentdepth-1]) ? true : false;
}

static void spse_nextLimdfsstate(const void *dfsconstinfo,
                                 DECLAREPTRDFSSTATE(aliasoutcol),
                                 unsigned long currentdepth,
                                 Uchar currentchar,
                                 GT_UNUSED const DECLAREPTRDFSSTATE(aliasincol))
{
  const Matchtaskinfo *mti = (const Matchtaskinfo *) dfsconstinfo;
  Limdfsstate *outcol = (Limdfsstate *) aliasoutcol;
#ifndef NDEBUG
  const Limdfsstate *incol = (const Limdfsstate *) aliasincol;
#endif

  gt_assert(ISNOTSPECIAL(currentchar));
  gt_assert(currentdepth > 0);
  gt_assert(incol->pathmatches);

  outcol->pathmatches = setpathmatch(mti->seedbitvector,
                                     mti->pattern,
                                     currentdepth,
                                     currentchar);
}

static void spse_inplacenextLimdfsstate(const void *dfsconstinfo,
                                        DECLAREPTRDFSSTATE(aliascol),
                                        unsigned long currentdepth,
                                        Uchar currentchar)
{
  const Matchtaskinfo *mti = (const Matchtaskinfo *) dfsconstinfo;
  Limdfsstate *col = (Limdfsstate *) aliascol;

  gt_assert(ISNOTSPECIAL(currentchar));
  gt_assert(currentdepth > 0);
  col->pathmatches = setpathmatch(mti->seedbitvector,
                                  mti->pattern,
                                  currentdepth,
                                  currentchar);
}

const AbstractDfstransformer *spse_AbstractDfstransformer(void)
{
  static const AbstractDfstransformer spse_adfst =
  {
    sizeof (Limdfsstate),
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
    NULL,
#ifdef SKDEBUG
    spse_showLimdfsstate,
#endif
  };
  return &spse_adfst;
}
