/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include <stdbool.h>
#include <ctype.h>
#include "libgtcore/env.h"
#include "libgtcore/hashtable.h"
#include "spacedef.h"
#include "dist-if.h"

 struct _Distribution
{
  Hashtable *hashdist;
};

static void freevalue(void *ptr,Env *env)
{
  env_error_check(env);
  FREESPACE(ptr);
}

Distribution *initdistribution(Env *env)
{
  Distribution *dist;

  env_error_check(env);
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

void addmultidistribution(Distribution *dist,unsigned long ind,
                          unsigned long howmany,Env *env)
{
  void *result;

  env_error_check(env);
  result = hashtable_get(dist->hashdist,(void *) ind);
  if (result == NULL)
  {
    unsigned long *newvalueptr;

    ALLOCASSIGNSPACE(newvalueptr,NULL,unsigned long,1);
    *newvalueptr = howmany;
    hashtable_add(dist->hashdist,(void *) ind,newvalueptr,env);
  } else
  {
    unsigned long *valueptr = (unsigned long *) result;

    (*valueptr) += howmany;
  }
}

void adddistribution(Distribution *dist,unsigned long ind,Env *env)
{
  env_error_check(env);
  addmultidistribution(dist,ind,(unsigned long) 1,env);
}

int foreachdistributionvalue(Distribution *dist,
                             int (*hashiter)(void *key, void *value,
                                             void *data, Env*),
                             void *data,Env *env)
{
  env_error_check(env);
  return hashtable_foreach(dist->hashdist,hashiter,data,env);
}
