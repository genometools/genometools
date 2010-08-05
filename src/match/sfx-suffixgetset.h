/*
  Copyright (c) 2010 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#ifndef SFX_SUFFIXGETSET_H
#define SFX_SUFFIXGETSET_H

#include "suffixptr.h"

typedef struct Suffixsortspace Suffixsortspace;

typedef void (*Dc_processunsortedrange)(void *,
                                        Suffixptr *,
                                        unsigned long,
                                        unsigned long,
                                        unsigned long);

Suffixsortspace *suffixsortspace_new(unsigned long numofentries);

Suffixsortspace *suffixsortspace_new_fromfile(int filedesc,
                                              const char *filename,
                                              unsigned long numofentries);

void suffixsortspace_delete(Suffixsortspace *suffixsortspace);

unsigned long suffixptrget(const Suffixsortspace *sssp,
                           const Suffixptr *subbucket,
                           unsigned long subbucketleft,
                           unsigned long idx);

void suffixptrset(Suffixsortspace *sssp,
                  Suffixptr *subbucket,
                  unsigned long subbucketleft,
                  unsigned long idx,
                  unsigned long value);

void suffixptrset2(const Suffixsortspace *sssp,
                   unsigned long idx,
                   unsigned long value);

unsigned long suffixptrget3(const Suffixsortspace *sssp,
                            unsigned long idx);

void suffixptrset3(Suffixsortspace *sssp,
                   unsigned long idx,
                   unsigned long value);

unsigned long gt_suffixsortspace_bucketleftidx_get(const Suffixsortspace *sssp);

void gt_suffixsortspace_bucketleftidx_set(Suffixsortspace *sssp,
                                          unsigned long value);

void gt_suffixsortspace_sortspace_set(Suffixsortspace *sssp,
                                      Suffixptr *sortspace);

unsigned long gt_suffixsortspace_offset_get(const Suffixsortspace *sssp);

void gt_suffixsortspace_sortspace_delete(Suffixsortspace *sssp);

void gt_suffixsortspace_offset_set(Suffixsortspace *sssp,
                                   unsigned long offset);

Suffixptr *gt_suffixsortspace_sortspace_get(const Suffixsortspace *sssp);

Suffixptr *gt_suffixsortspace_leftadjust(const Suffixsortspace *sssp);

#endif
