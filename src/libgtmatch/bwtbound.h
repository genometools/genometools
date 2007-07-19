#ifndef BWTBOUND_H
#define BWTBOUND_H

#include "types.h"
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
