/*
  Copyright (c) 2006 Gordon Gremme <gremme@@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FPTR_H
#define FPTR_H

/* the generic function pointers */
typedef int  (*Fptr)();
typedef int  (*Compare)(const void*, const void*);
typedef void (*Free)(void*);

#endif
