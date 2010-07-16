/*
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#include "example_rep.h"       /* we need access to the class struct */
#include "core/ma.h"           /* we need to allocate memory */

GtExample* gt_example_create(const GtExampleClass *ec)
{
  GtExample *e = gt_calloc(1, ec->size);   /* allocate memory */
  e->c_class = ec;                         /* assign interface */
  return e;
}

int gt_example_run(GtExample *e)
{
  gt_assert(e && e->c_class && e->c_class->run);
  return e->c_class->run(e);               /* call implementation-specific
                                              function */
}

void gt_example_delete(GtExample *e)
{
  if (!e) return;
  gt_assert(e && e->c_class);
  if (e->c_class->delete != NULL) {
    e->c_class->delete(e);                 /* delete implementation-specific
                                              members */
  }
  gt_free(e);                              /* delete interface */
}
