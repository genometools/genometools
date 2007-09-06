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

#ifndef MAPSPECDEF_H
#define MAPSPECDEF_H

#include "libgtcore/env.h"
#include "arraydef.h"

#define NEWMAPSPEC(PTR,TYPE,ELEMS)\
        GETNEXTFREEINARRAY(mapspecptr,mapspectable,Mapspecification,10);\
        mapspecptr->typespec = TYPE ## Type;\
        mapspecptr->startptr = &(PTR);\
        mapspecptr->sizeofunit = sizeof (TYPE);\
        mapspecptr->numofunits = ELEMS;\
        mapspecptr->name = #PTR

typedef enum
{
  UcharType,
  UshortType,
  Uint32Type,
  Uint64Type,
  BitstringType,
  SeqposType,
  BwtboundType,
  PairBwtidxType
} Typespec;

typedef struct
{
  Typespec typespec;
  char *name;
  void *startptr;
  size_t sizeofunit;
  unsigned long numofunits;
} Mapspecification;

DECLAREARRAYSTRUCT(Mapspecification);

typedef void(*Assignmapspec)(ArrayMapspecification *,void *,Env *);

#endif
