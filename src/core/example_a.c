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
#include "example_a.h"
#include "example_rep.h"

struct GtExampleA {
  GtExample parent_instance;
  unsigned long my_property;
};

static int gt_example_a_run(GtExample *e) /* hidden from outside  */
{
  GtExampleA *ea = (GtExampleA*) e;       /* downcast to specific type */
  printf("%lu", ea->my_property);         /* run functionality */
  return 0;
}

/* map static local method to interface */
const GtExampleClass* gt_example_a_class(void)
{
  static const GtExampleClass ec = { sizeof (GtExampleA),
                                     gt_example_a_run,
                                     NULL };
  return &ec;
}

GtExample* gt_example_a_new(void)
{
  GtExample *e = gt_example_create(gt_example_a_class());
  GtExampleA *ea = (GtExampleA*) e;       /* downcast to specific type */
  ea->my_property = 3;                    /* access private implementation
                                             member */
  return e;
}
