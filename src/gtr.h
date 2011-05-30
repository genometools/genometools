/*
  Copyright (c) 2003-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef GTR_H
#define GTR_H

#include <stdio.h>
#include "core/init_api.h"
#include "core/error.h"
#include "core/option_api.h"

/* The GenomeTools runtime (gtr) */
typedef struct GtR GtR;

GtR*     gtr_new(GtError*);
GtOPrval gtr_parse(GtR*, int *parsed_args, int argc, const char **argv,
                   GtError*);
void     gtr_register_components(GtR*);
int      gtr_run(GtR*, int argc, const char **argv, GtError*);
void     gtr_delete(GtR*);

#endif
