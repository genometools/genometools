/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdbool.h>
#include <ctype.h>
#include "libgtcore/env.h"
#include "libgtcore/hashtable.h"
#include "types.h"
#include "spacedef.h"
#include "dist-if.h"

 struct _Distribution
{
  Hashtable *hashdist;
};

static void freevalue(void *ptr,Env *env)
{
  FREESPACE(ptr);
}

Distribution *initdistribution(Env *env)
{
  Distribution *dist;
  ALLOCASSIGNSPACE(dist,NULL,Distribution,1);
  dist->hashdist = hashtable_new(HASH_DIRECT, NULL, freevalue, env);
  return dist;
}

void freedistribution(Distribution **dist,Env *env)
{
  hashtable_delete((*dist)->hashdist,env);
  FREESPACE(*dist);
}

/* XXX: allow howmany to be of type Uint64 */

void addmultidistribution(Distribution *dist,Uint ind,Uint howmany,Env *env)
{
  void *result;

  result = hashtable_get(dist->hashdist,(void *) ind);
  if(result == NULL)
  {
    Uint *newvalueptr;

    ALLOCASSIGNSPACE(newvalueptr,NULL,Uint,1);
    *newvalueptr = howmany;
    hashtable_add(dist->hashdist,(void *) ind,newvalueptr,env);
  } else
  {
    Uint *valueptr = (Uint *) result;

    (*valueptr) += howmany;
  }
}

void adddistribution(Distribution *dist,Uint ind,Env *env)
{
  addmultidistribution(dist,ind,UintConst(1),env);
}

int foreachdistributionvalue(Distribution *dist,
                             int (*hashiter)(void *key, void *value,
                                             void *data, Env*),
                             void *data,Env *env)
{
  return hashtable_foreach(dist->hashdist,hashiter,data,env);
}
