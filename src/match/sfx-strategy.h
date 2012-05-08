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

#include <inttypes.h>
#include <stdbool.h>

#define MAXINSERTIONSORTDEFAULT  3UL
#define MAXBLTRIESORTDEFAULT     1000UL
#define MAXCOUNTINGSORTDEFAULT   4000UL

typedef struct
{
  unsigned long maxwidthrealmedian,
                maxinsertionsort,
                maxbltriesort,
                maxcountingsort;
  unsigned int differencecover,
               userdefinedsortmaxdepth,
               spmopt_minlength;
  bool cmpcharbychar, /* compare suffixes character by character instead
                         of comparing entire words (only for two bit
                         encoding) */
       storespecialcodes,
       iteratorbasedkmerscanning,
       suftabuint,
       onlybucketinsertion,
       kmerswithencseqreader,
       dccheck,
       samplewithprefixlengthnull,
       noshortreadsort,
       outsuftabonfile,
       withradixsort;
} Sfxstrategy;

 /*@unused@*/ static inline void defaultsfxstrategy(Sfxstrategy *sfxstrategy,
                                                    bool cmpcharbychar)
{
  sfxstrategy->maxwidthrealmedian = 1UL;
  sfxstrategy->maxinsertionsort = MAXINSERTIONSORTDEFAULT;
  sfxstrategy->maxbltriesort = MAXBLTRIESORTDEFAULT;
  sfxstrategy->maxcountingsort = MAXCOUNTINGSORTDEFAULT;
  sfxstrategy->differencecover = 0;
  sfxstrategy->cmpcharbychar = cmpcharbychar;
  sfxstrategy->spmopt_minlength = 0;
  sfxstrategy->storespecialcodes = false;
  sfxstrategy->iteratorbasedkmerscanning = false;
  sfxstrategy->suftabuint = false;
  sfxstrategy->onlybucketinsertion = false;
  sfxstrategy->kmerswithencseqreader = false;
  sfxstrategy->dccheck = false;
  sfxstrategy->samplewithprefixlengthnull = false;
  sfxstrategy->outsuftabonfile = true;
  sfxstrategy->noshortreadsort = false;
  sfxstrategy->withradixsort = false;
  sfxstrategy->userdefinedsortmaxdepth = 0;
}

#endif
