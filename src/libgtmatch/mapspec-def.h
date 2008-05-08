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

#ifndef MAPSPEC_DEF_H
#define MAPSPEC_DEF_H

#include "libgtcore/arraydef.h"
#include "libgtcore/str.h"
#include "libgtcore/error.h"

#define NEWMAPSPEC(PTR,TYPE,ELEMS)\
        GETNEXTFREEINARRAY(mapspecptr,mapspectable,Mapspecification,10);\
        mapspecptr->typespec = TYPE ## Type;\
        mapspecptr->startptr = &(PTR);\
        mapspecptr->sizeofunit = sizeof (TYPE);\
        mapspecptr->numofunits = ELEMS;\
        mapspecptr->name = #PTR

typedef unsigned long Unsignedlong;

typedef enum
{
  UcharType,
  UshortType,
  Uint32Type,
  Uint64Type,
  UnsignedlongType,
  BitstringType,
  SeqposType,
  BwtboundType,
  PairBwtidxType,
  TwobitencodingType
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

typedef void(*Assignmapspec)(ArrayMapspecification *,void *,bool);

int fillmapspecstartptr(Assignmapspec assignmapspec,
                        void **mappeduserptr,
                        void *assignmapinfo,
                        const Str *tmpfilename,
                        unsigned long expectedsize,
                        Error *err);

int flushtheindex2file(FILE *fp,
                       Assignmapspec assignmapspec,
                       void *assignmapinfo,
                       unsigned long expectedsize,
                       Error *err);

#endif
