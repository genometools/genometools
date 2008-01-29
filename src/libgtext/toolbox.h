/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#include "libgtcore/error.h"
#include "libgtcore/tool.h"

/* the toolbox class */
typedef struct Toolbox Toolbox;

typedef int (*Toolfunc)(int argc, const char **argv, Error*);

Toolbox* toolbox_new(void);

/* Add <tool> with name <toolname> to <toolbox>. Takes ownership of <tool>. */
void     toolbox_add_tool(Toolbox*, const char *toolname, Tool *tool);
/* Get Tool with name <toolname> from <toolbox>. */
Tool*    toolbox_get_tool(Toolbox*, const char *toolname);
bool     toolbox_has_tool(const Toolbox*, const char *toolname); /*deprecated */
void     toolbox_add(Toolbox*, const char *toolname, Toolfunc); /* deprecated */
Toolfunc toolbox_get(const Toolbox*, const char *toolname); /* deprecated */
/* shows all tools except tools with toolname ``dev'' */
int      toolbox_show(/*@unused@*/ const char *progname, void *toolbox, Error*);
void     toolbox_delete(Toolbox*);

#endif
