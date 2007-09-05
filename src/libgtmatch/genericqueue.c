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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "libgtcore/env.h"
#include "spacedef.h"
#include "gqueue-def.h"

typedef struct _Genericelementofqueue
{
  void *contents;
  struct _Genericelementofqueue *previous, *next;
} Genericelementofqueue;

 struct Genericqueue
{
  Genericelementofqueue *head,
                        *tail;
  unsigned long numberofelements;
};

/*
  This file defines a generic type for queues.
*/

Genericqueue *emptyqueuegeneric(Env *env)
{
  Genericqueue *q;

  env_error_check(env);
  ALLOCASSIGNSPACE(q,NULL,Genericqueue,1);
  assert(q != NULL);
  q->numberofelements = 0;
  q->head = NULL;
  q->tail = NULL;
  assert(q != NULL);
  return q;
}

unsigned long sizeofgenericqueue(const Genericqueue *q)
{
  return q->numberofelements;
}

bool queueisemptygeneric(const Genericqueue *q)
{
  if (q->numberofelements == 0)
  {
    return true;
  }
  return false;
}

void enqueuegeneric(Genericqueue *q,void *contents,Env *env)
{
  Genericelementofqueue *newqueueelem;

  env_error_check(env);
  ALLOCASSIGNSPACE(newqueueelem,NULL,Genericelementofqueue,1);
  newqueueelem->contents = contents;
  newqueueelem->previous = NULL;
  newqueueelem->next = q->tail;
  if (q->numberofelements == 0)
  {
    q->head = newqueueelem;
  } else
  {
    q->tail->previous = newqueueelem;
  }
  q->tail = newqueueelem;
  q->numberofelements++;
}

/*@null@*/ void *dequeuegeneric(Genericqueue *q,Env *env)
{
  void *contents;
  Genericelementofqueue *oldheadptr;

  env_error_check(env);
  if (q->numberofelements == 0)
  {
    env_error_set(env,"dequeuegeneric(emptyqueue) is undefined");
    return NULL;
  }
  oldheadptr = q->head;
  q->head = q->head->previous;
  if (q->head == NULL)
  {
    q->tail = NULL;
  } else
  {
    q->head->next = NULL;
  }
  contents = oldheadptr->contents;
  FREESPACE(oldheadptr);
  q->numberofelements--;
  return contents;
}

/*@null@*/ void *headofqueuegeneric(const Genericqueue *q)
{
  if (q->numberofelements == 0)
  {
    return NULL;
  }
  return q->head->contents;
}

int overallqueuelementsgeneric(Genericqueue *q,
                               GenericQueueprocessor queueprocessor,
                               void *info)
{
  Genericelementofqueue *current;

  if (q->numberofelements > 0)
  {
    for (current = q->head;
        current != NULL;
        current = current->previous)
    {
      if (queueprocessor(current->contents,info) != 0)
      {
        return -1;
      }
    }
  }
  return 0;
}

void wrapqueuegeneric(bool freecontents,Genericqueue **q,Env *env)
{
  while (!queueisemptygeneric(*q))
  {
    if (freecontents)
    {
      FREESPACE((*q)->head->contents);
    }
    (void) dequeuegeneric(*q,env);
  }
  FREESPACE(*q);
}
