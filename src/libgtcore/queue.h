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

#ifndef QUEUE_H
#define QUEUE_H

#include <stdio.h>
#include "libgtcore/error.h"

typedef struct Queue Queue;

typedef int (*QueueProcessor)(void *elem, void *info, Error*);

Queue*        queue_new(void);
void          queue_add(Queue*, void*);
void*         queue_get(Queue*);
void*         queue_head(Queue*);
int           queue_iterate(Queue*, QueueProcessor, void *info, Error*);
unsigned long queue_size(const Queue*);
void          queue_delete_with_contents(Queue*);
void          queue_delete(Queue*);

#endif
