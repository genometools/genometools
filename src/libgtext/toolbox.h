/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef TOOLBOX_H
#define TOOLBOX_H

#include <gtcore.h>

/* the toolbox class */
typedef struct Toolbox Toolbox;

typedef int (*Tool)(int argc, const char **argv, Env*);

Toolbox* toolbox_new(Env*);
void     toolbox_add(Toolbox*, const char *toolname, Tool, Env*);
Tool     toolbox_get(const Toolbox*, const char *toolname);
int      toolbox_show(/*@unused@*/ const char *progname, void *toolbox, Env*);
void     toolbox_delete(Toolbox*, Env*);

#endif
