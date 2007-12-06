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

#ifndef TEST_PAIRWISE_H
#define TEST_PAIRWISE_H

#include <stdbool.h>
#include "libgtcore/symboldef.h"

typedef void (*Checkcmppairfuntype)(bool,
                                    const Uchar *,unsigned long,
                                    const Uchar *,unsigned long);

void runcheckfunctionontwofiles(Checkcmppairfuntype checkfunction,
                               const char *file1,
                               const char *file2);

unsigned long runcheckfunctionontext(Checkcmppairfuntype checkfunction,
                                     const char *text);

unsigned long applycheckfunctiontotext(const Uchar *text,
                                       unsigned long textlen,
                                       void *info);

unsigned long runcheckfunctiononalphalen(Checkcmppairfuntype checkfunction,
                                         const char *charlist,
                                         unsigned long len);

void checkgreedyunitedist(/*@unused@*/ bool forward,
                          const Uchar *useq,
                          unsigned long ulen,
                          const Uchar *vseq,
                          unsigned long vlen);

#endif
