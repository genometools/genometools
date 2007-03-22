/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinforfatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FA_H
#define FA_H

#include <stdio.h>
#include <stdlib.h>

/* the file allocator class */
typedef struct FA FA;

FA*     fa_new(Env*);

/* check if all allocated file pointer have been released, prints to stderr */
int     fa_check_fptr_leak(FA*, Env*);
/* check if all allocated memory maps have been freed, prints to stderr */
int     fa_check_mmap_leak(FA*, Env*);
void    fa_delete(FA*, Env*);

#endif
