/*
  Copyright (c) 2003-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GTR_H
#define GTR_H

#include <stdio.h>
#include "libgtcore/env.h"
#include "libgtcore/option.h"

/* The GenomeTools runtime (gtr) */
typedef struct GTR GTR;

GTR*   gtr_new(Env*);
OPrval gtr_parse(GTR*, int *parsed_args, int argc, const char **argv, Env*);
void   gtr_register_components(GTR*, Env*);
int    gtr_run(GTR*, int argc, const char **argv, Env*);
void   gtr_delete(GTR*, Env*);

#endif
