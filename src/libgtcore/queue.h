/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef QUEUE_H
#define QUEUE_H

#include <stdio.h>
#include "libgtcore/env.h"

/* A simple queue implementation (based on arrays). Only memory efficient if
   the queue becomes empty regularly during usage (see implementation for
   details).  */
typedef struct Queue Queue;

Queue*        queue_new(size_t, Env*);
void*         queue_get(Queue*);
void*         queue_get_elem(Queue*, unsigned long);
#define       queue_add(q, elem, env)\
              queue_add_elem(q, &(elem), sizeof (elem), env)
void          queue_add_elem(Queue*, void*, size_t, Env*);
unsigned long queue_size(const Queue*);
void          queue_delete(Queue*, Env*);

#endif
