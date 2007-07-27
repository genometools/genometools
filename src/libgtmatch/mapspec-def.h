/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
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
