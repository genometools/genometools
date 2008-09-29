/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include "core/array.h"
#include "core/class_alloc.h"
#include "core/ma.h"

static GtArray *c_classes = NULL;

void* gt_class_alloc(size_t size)
{
  void *c_class;
  if (!c_classes)
    c_classes = gt_array_new(sizeof (void*));
  c_class = gt_calloc(1, size);
  gt_array_add(c_classes, c_class);
  return c_class;
}

void  gt_class_alloc_clean(void)
{
  unsigned long i;
  if (!c_classes) return;
  for (i = 0; i < gt_array_size(c_classes); i++)
    gt_free(*(void**) gt_array_get(c_classes, i));
  gt_array_delete(c_classes);
  c_classes = NULL;
}
