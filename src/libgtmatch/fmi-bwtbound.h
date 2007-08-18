/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef BWTBOUND_H
#define BWTBOUND_H

#include "seqpos-def.h"
#include "arraydef.h"

typedef struct
{
  Seqpos lbound,
         ubound;
} Bwtbound;

typedef struct
{
  Seqpos bwtpos,
         suftabvalue;
} PairBwtidx;

DECLAREARRAYSTRUCT(PairBwtidx);

#endif
