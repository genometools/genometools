/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef TOOLBOX_H
#define TOOLBOX_H

#include "libgtcore/env.h"

/* the toolbox class */
typedef struct Toolbox Toolbox;

typedef int (*Tool)(int argc, const char **argv, Env*);

Toolbox* toolbox_new(void);
void     toolbox_add(Toolbox*, const char *toolname, Tool);
Tool     toolbox_get(const Toolbox*, const char *toolname);
/* shows all tools except tools with toolname ``dev'' */
int      toolbox_show(/*@unused@*/ const char *progname, void *toolbox, Error*);
void     toolbox_delete(Toolbox*);

#endif
