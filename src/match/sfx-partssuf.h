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

#ifndef SFX_PARTSSUF_H
#define SFX_PARTSSUF_H

#include "core/logger.h"
#include "core/codetype.h"
#include "bcktab.h"

typedef struct Suftabparts Suftabparts;

Suftabparts *gt_newsuftabparts(unsigned int numofparts,
                               const GtBcktab *bcktab,
                               unsigned long numofsuffixestoinsert,
                               unsigned long fullspecials,
                               GtLogger *logger);

GtCodetype stpgetcurrentmincode(unsigned int part,
                              const Suftabparts *suftabparts);

GtCodetype stpgetcurrentmaxcode(unsigned int part,
                                const Suftabparts *suftabparts);

unsigned long stpgetcurrentsuftaboffset(unsigned int part,
                                        const Suftabparts *suftabparts);

unsigned long stpgetcurrentsumofwdith(unsigned int part,
                               const Suftabparts *suftabparts);

unsigned long stpgetcurrentwidthofpart(unsigned int part,
                                const Suftabparts *suftabparts);

unsigned long stpgetlargestsuftabwidth(const Suftabparts *suftabparts);

unsigned int stpgetnumofparts(const Suftabparts *suftabparts);

unsigned long stpgetlargestsizemappedpartwise(const Suftabparts *suftabparts);

unsigned long stpnumofsuffixestoinsert(const Suftabparts *suftabparts);

void gt_freesuftabparts(Suftabparts *suftabparts);

#endif
