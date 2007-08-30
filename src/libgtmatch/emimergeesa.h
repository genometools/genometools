/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef EMIMERGEESA_H
#define EMIMERGEESA_H

#include "seqpos-def.h"
#include "sarr-def.h"
#include "trieins-def.h"

#define SIZEOFMERGERESULTBUFFER BUFSIZ

typedef struct
{
  unsigned int idx;  /* index of genome in list of all genomes */
  Seqpos startpos;   /* in the range [0..totallength single index] */
} Indexedsuffix;

typedef struct
{
  unsigned int nextaccessidx,  /* in the range [0..SIZEOFMERGERESULTBUFFER] */
               nextstoreidx;   /* in the range [0..SIZEOFMERGERESULTBUFFER] */
  Seqpos lcptabstore[SIZEOFMERGERESULTBUFFER];
  Indexedsuffix suftabstore[SIZEOFMERGERESULTBUFFER];
  bool lastpage;
} Suflcpbuffer;

typedef struct
{
  uint64_t ident;              /* can be arbitrary large */
  unsigned int numofentries,   /* in the range [0..numofindexes-1] */
               numofindexes;   /* number of indexes */
  Seqpos *nextpostable;        /* in the range [0..totallength single index] */
  Suflcpbuffer buf;
  Trierep trierep;
  Suffixarray *suffixarraytable;
  Alphabet *alpha;
} Emissionmergedesa;

#endif
