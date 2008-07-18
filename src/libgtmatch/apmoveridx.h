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

#ifndef APMOVERIDX_H
#define APMOVERIDX_H
#include "libgtcore/unused.h"
#include "defined-types.h"

typedef unsigned long Aliasdfsstate;

#define DECLAREPTRDFSSTATE(V)\
        Aliasdfsstate * V

#undef SKDEBUG

typedef struct
{
  size_t sizeofdfsstate;
  void (*initdfsconstinfo)(void *dfsconstinfo,
                           unsigned int alphasize,
                           const Uchar *pattern,
                           unsigned long patternlength,
                           unsigned long maxdistance,
                           Seqpos maxintervalwidth);
  void *(*allocatedfsconstinfo)(unsigned int alphasize);
  void (*freedfsconstinfo)(void **dfsconstinfo);
  unsigned long (*limdfsnextstep)(const DECLAREPTRDFSSTATE(aliascolumn),
                                  Seqpos width,
                                  const void *dfsconstinfo);
  void (*nextDfsstate)(const void *dfsconstinfo,
                       DECLAREPTRDFSSTATE(aliasoutstate),
                       UNUSED unsigned long startscore,
                       Uchar currentchar,
                       const DECLAREPTRDFSSTATE(aliasinstate));
  void (*inplacenextDfsstate)(const void *dfsconstinfo,
                              DECLAREPTRDFSSTATE(aliasstate),
                              Uchar currentchar);
  void (*initLimdfsstate)(DECLAREPTRDFSSTATE(aliasstate),
                          const void *dfsconstinfo);
#ifdef SKDEBUG
  void (*showLimdfsstate)(const DECLAREPTRDFSSTATE(aliasstate),
                          unsigned long score,
                          const void *dfsconstinfo);
#endif
} AbstractDfstransformer;

const AbstractDfstransformer *apm_AbstractDfstransformer(void);

Definedunsignedlong apm_findshortestmatch(const Encodedsequence *encseq,
                                          bool nospecials,
                                          unsigned int alphasize,
                                          const Uchar *pattern,
                                          unsigned long patternlength,
                                          unsigned long maxdistance,
                                          Seqpos startpos);

#endif
