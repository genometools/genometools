/*
  Copyright (c) 2003-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GTR_H
#define GTR_H

#include <stdio.h>

/* The GenomeTools runtime (gtr) */
typedef struct GTR GTR;

GTR* gtr_new(void);
int  gtr_parse(GTR*, int argc, char **argv);
void gtr_register_components(GTR*);
void gtr_show_help(GTR*, FILE*);
int  gtr_run(GTR*, int argc, char **argv);
void gtr_free(GTR*);

#endif
