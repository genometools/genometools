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

#ifndef SFX_STRATEGY_H
#define SFX_STRATEGY_H

#include <stdbool.h>
#include "core/defined-types.h"

#define MAXINSERTIONSORTDEFAULT  3UL
#define MAXBLTRIESORTDEFAULT     1000UL
#define MAXCOUNTINGSORTDEFAULT   4000UL

typedef struct
{
  Definedunsignedint ssortmaxdepth;
  unsigned long maxwidthrealmedian,
                maxcountingsort,
                maxinsertionsort,
                maxbltriesort;
  unsigned int differencecover;
  bool cmpcharbychar, /* compare suffixes character by character instead
                         of comparing entire words (only for two bit
                         encoding) */
       storespecialcodes,
       streamsuftab,
       absoluteinversesuftab,
       hashexceptions,
       iteratorbasedkmerscanning;
} Sfxstrategy;

 /*@unused@*/ static inline void defaultsfxstrategy(Sfxstrategy *sfxstrategy,
                                                    bool cmpcharbychar)
{
  sfxstrategy->ssortmaxdepth.defined = false;
  sfxstrategy->maxwidthrealmedian = 1UL;
  sfxstrategy->maxcountingsort = MAXCOUNTINGSORTDEFAULT;
  sfxstrategy->maxinsertionsort = MAXINSERTIONSORTDEFAULT;
  sfxstrategy->maxbltriesort = MAXBLTRIESORTDEFAULT;
  sfxstrategy->differencecover = 0;
  sfxstrategy->cmpcharbychar = cmpcharbychar;
  sfxstrategy->storespecialcodes = false;
  sfxstrategy->streamsuftab = false;
  sfxstrategy->absoluteinversesuftab = false;
  sfxstrategy->hashexceptions = false;
  sfxstrategy->iteratorbasedkmerscanning = false;
}

#endif
