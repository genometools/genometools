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

#ifndef ABSDFSTRANS_IMP_H
#define ABSDFSTRANS_IMP_H
#include "core/defined-types.h"
#include "absdfstrans-def.h"
#include "procmatch.h"

typedef unsigned long Aliasdfsstate;

#define DECLAREPTRDFSSTATE(V)\
        Aliasdfsstate * V

#undef SKDEBUG

typedef enum
{
  Limdfssuccess,  /* success of traversal */
  Limdfscontinue, /* no success, but still have the chance to find result */
  Limdfsstop      /* no success possible */
} Limdfsstatus;

typedef struct
{
  Limdfsstatus status;
  unsigned long pprefixlen, /* only defined if limdfsstatus = Limdfssuccess */
                distance;   /* only defined if limdfsstatus = Limdfssuccess */
} Limdfsresult;

struct AbstractDfstransformer
{
  size_t sizeofdfsstate;
  Limdfsconstinfo *(*allocatedfsconstinfo)(unsigned int alphasize);
  void (*initdfsconstinfo)(Limdfsconstinfo *dfsconstinfo,
                           unsigned int alphasize,
                           ...);
  void (*extractdfsconstinfo)(Processresult processresult,
                              void *processinfo,
                              const void *patterninfo,
                              Limdfsconstinfo *dfsconstinfo);
  void (*freedfsconstinfo)(Limdfsconstinfo **dfsconstinfo);
  void (*initrootLimdfsstate)(DECLAREPTRDFSSTATE(aliasstate),
                              Limdfsconstinfo *dfsconstinfo);
  void (*initLimdfsstackelem)(DECLAREPTRDFSSTATE(aliasstate));
  void (*freeLimdfsstackelem)(DECLAREPTRDFSSTATE(aliasstate));
  void (*copyLimdfsstate)(DECLAREPTRDFSSTATE(deststate),
                          const DECLAREPTRDFSSTATE(srcstate),
                          Limdfsconstinfo *dfsconstinfo);
  void (*fullmatchLimdfsstate)(Limdfsresult *limdfsresult,
                               DECLAREPTRDFSSTATE(aliascolumn),
                               unsigned long left,
                               unsigned long right,
                               unsigned long width,
                               unsigned long currentdepth,
                               Limdfsconstinfo *dfsconstinfo);
  void (*nextLimdfsstate)(const Limdfsconstinfo *dfsconstinfo,
                          DECLAREPTRDFSSTATE(aliasoutstate),
                          unsigned long currentdepth,
                          GtUchar currentchar,
                          const DECLAREPTRDFSSTATE(aliasinstate));
  void (*inplacenextLimdfsstate)(const Limdfsconstinfo *dfsconstinfo,
                                 DECLAREPTRDFSSTATE(aliasstate),
                                 unsigned long currentdepth,
                                 GtUchar currentchar);
#ifdef SKDEBUG
  void (*showLimdfsstate)(const DECLAREPTRDFSSTATE(aliasstate),
                          unsigned long currentdepth,
                          const Limdfsconstinfo *dfsconstinfo);
#endif
};

#endif
