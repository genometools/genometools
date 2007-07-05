/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdio.h>
#include <time.h>
#include "libgtcore/env.h"
#include "measure-time-if.h"

 struct _Measuretime
{
  clock_t startclock, overalltime;
  const char *eventdescription;
};

void inittheclock(Measuretime **mtime,const char *event,Env *env)
{
  *mtime = env_ma_malloc(env,sizeof (Measuretime));
  (*mtime)->startclock = clock();
  (*mtime)->overalltime = 0;
  (*mtime)->eventdescription = event;
}

void deliverthetime(FILE *fp,Measuretime *mtime,const char *newevent,Env *env)
{
  clock_t stopclock;

  stopclock = clock();
  fprintf(fp,"TIME %s %.2f\n",mtime->eventdescription,
             (double) (stopclock-mtime->startclock)/(double) CLOCKS_PER_SEC);
  (void) fflush(fp);
  mtime->overalltime += (stopclock - mtime->startclock);
  if (newevent == NULL)
  {
    fprintf(fp,"TIME overall %.2f\n",
                (double) mtime->overalltime/(double) CLOCKS_PER_SEC);
    (void) fflush(fp);
    env_ma_free(mtime,env);
  } else
  {
    mtime->startclock = stopclock;
    mtime->eventdescription = newevent;
  }
}
