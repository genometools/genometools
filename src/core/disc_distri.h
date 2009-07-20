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

#ifndef DISC_DISTRI_H
#define DISC_DISTRI_H

#include "core/error.h"
#include "core/file.h"

/* A discrete distribution */
typedef struct GtDiscDistri GtDiscDistri;

typedef void (*GtDiscDistriIterFunc)(unsigned long key,
                                     unsigned long long value,
                                     void *data);

GtDiscDistri*        gt_disc_distri_new(void);
void               gt_disc_distri_add(GtDiscDistri*, unsigned long);
void               gt_disc_distri_add_multi(GtDiscDistri*, unsigned long,
                                        unsigned long long);
unsigned long long gt_disc_distri_get(const GtDiscDistri*, unsigned long);
void               gt_disc_distri_show(const GtDiscDistri*); /* on stdout */
void               gt_disc_distri_show_generic(const GtDiscDistri*, GtGenFile*);
void               gt_disc_distri_foreach(const GtDiscDistri*,
                                          GtDiscDistriIterFunc,
                                          void *data);
int                gt_disc_distri_unit_test(GtError*);
void               gt_disc_distri_delete(GtDiscDistri*);

#endif
