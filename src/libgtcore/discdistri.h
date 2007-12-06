/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c)      2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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

#ifndef DISCDISTRI_H
#define DISCDISTRI_H

#include "libgtcore/env.h"
#include "libgtcore/genfile.h"

/* A discrete distribution */
typedef struct DiscDistri DiscDistri;

typedef void (*DiscDistriIterFunc)(unsigned long key, unsigned long long value,
                                   void *data);

DiscDistri*        discdistri_new(void);
void               discdistri_add(DiscDistri*, unsigned long);
void               discdistri_add_multi(DiscDistri*, unsigned long,
                                        unsigned long long);
unsigned long long discdistri_get(const DiscDistri*, unsigned long);
void               discdistri_show(const DiscDistri*); /* on stdout */
void               discdistri_show_generic(const DiscDistri*, GenFile*);
void               discdistri_foreach(const DiscDistri*, DiscDistriIterFunc,
                                      void *data);
int                discdistri_unit_test(Error*);
void               discdistri_delete(DiscDistri*);

#endif
