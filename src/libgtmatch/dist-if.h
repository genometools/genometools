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

#ifndef DIST_IF_H
#define DIST_IF_H
#include "libgtcore/env.h"

typedef struct Distribution Distribution;

Distribution *initdistribution(Env *env);

void freedistribution(Distribution **dist,Env *env);

void addmultidistribution(Distribution *dist,unsigned long ind,
                          unsigned long howmany,Env *env);

void adddistribution(Distribution *dist,unsigned long ind,Env *env);

int foreachdistributionvalue(Distribution *dist,
                             int (*hashiter)(void *key, void *value,
                                             void *data, Env*),
                             void *data,Env *env);

#endif
