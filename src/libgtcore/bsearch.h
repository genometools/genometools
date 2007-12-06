/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef BSEARCH_H
#define BSEARCH_H

#include "libgtcore/array.h"
#include "libgtcore/bittab.h"

typedef int (*Compar)(const void*, const void*, const void *data);

/* similar interface to bsearch(3), except that the Compar function gets an
   additional <data> pointer */
void* bsearch_data(const void *key, const void *base, size_t nmemb, size_t size,
                   Compar, void *data);

/* similar interface to bsearch_data(), except that all members which compare as
   equal are stored in the <members> array. The order in which the elements
   are added is undefined */
void  bsearch_all(Array *members, const void *key, const void *base,
                  size_t nmemb, size_t size, Compar, void *data);

/* similar interface to bsearch_all(). Additionally, if a bittab is given (which
   must be of size <nmemb>), the bits corresponding to the found elements are
   marked (i.e., set) */
void  bsearch_all_mark(Array *members, const void *key, const void *base,
                       size_t nmemb, size_t size, Compar, void *data, Bittab*);

int   bsearch_unit_test(Error*);

#endif
