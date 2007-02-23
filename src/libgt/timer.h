/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef TIMER_H
#define TIMER_H

#include <sys/time.h>
#include <sys/resource.h>
#include <stdio.h>
#include "env.h"

/* the timer class */
typedef struct Timer Timer;

Timer* timer_new(Env*);
void   timer_start(Timer*);
void   timer_stop(Timer*);
void   timer_show(Timer*, FILE*);
void   timer_del(Timer*, Env*);

#endif
