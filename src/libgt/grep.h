/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GREP_H
#define GREP_H

#include <sys/types.h>
#include <assert.h>
#include <regex.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <libgt/env.h>

/* sets 'match' to true if pattern matches line, to false otherwise */
int  grep(bool *match, const char *pattern, const char *line, Env*);
int  grep_unit_test(Env*);

#endif
