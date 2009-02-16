/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#ifndef IDXLOCALIDP_H
#define IDXLOCALIDP_H
#include "core/symboldef.h"
#include "seqpos-def.h"
#include "absdfstrans-def.h"

const AbstractDfstransformer *locali_AbstractDfstransformer(void);

typedef struct Localitracebackstate Localitracebackstate;

Localitracebackstate *newLocalitracebackstate(const Uchar *characters,
                                              Uchar wildcardshow);

void reinitLocalitracebackstate(Localitracebackstate *tbs,
                                Seqpos dbprefixlen,unsigned long pprefixlen);

void processelemLocalitracebackstate(Localitracebackstate *tbs,
                                     Uchar currentchar,
                                     const void *aliasstate);

void showLocalitracebackstate(const void *dfsconstinfo,
                              const Localitracebackstate *tbs);

void freeLocalitracebackstate(Localitracebackstate *);

#endif
