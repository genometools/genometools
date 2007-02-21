/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GTDATA_H
#define GTDATA_H

#include "str.h"

/* This module defines functions working on the ``gtdata'' directory */

/* get the path to the gtdata/ directory (including it) for the given 'prog' */
Str* gtdata_get_path(const char *prog, Env*);

/* execute helpfile gtdata/doc/progname.lua */
int gtdata_show_help(const char *progname, void *unused, Env*);

#endif
