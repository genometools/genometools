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

#include <stdio.h>
#include "example_b.h"
#include "example_rep.h"
#include "core/str.h"

struct GtExampleB {
  GtExample parent_instance;
  GtStr *my_property;       /* this is an allocated object we need to
                               memory-manage */
};

static int gt_example_b_run(GtExample *e) /* hidden from outside  */
{
  GtExampleB *eb = (GtExampleB*) e;          /* downcast to specific type */
  printf("%s", gt_str_get(eb->my_property)); /* run functionality */
  return 0;
}

static void gt_example_b_delete(GtExample *e)
{
  GtExampleB *eb = (GtExampleB*) e;       /* downcast to specific type */
  gt_str_delete(eb->my_property);
}

/* map static local method to interface */
const GtExampleClass* gt_example_b_class(void)
{
  static const GtExampleClass ec = { sizeof (GtExampleB),
                                     gt_example_b_run,
                                     gt_example_b_delete };
  return &ec;
}

GtExample* gt_example_b_new(void)
{
  GtExample *e = gt_example_create(gt_example_b_class());
  GtExampleB *eb = (GtExampleB*) e;       /* downcast to specific type */
  /* initialize private implementation member */
  eb->my_property = gt_str_new_cstr("This is Example B!");
  return e;
}
