/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef SARR_DEF_H
#define SARR_DEF_H

#include "types.h"
#include "alphadef.h"
#include "encseq-def.h"

typedef struct
{
  uint32_t numofdbsequences;
  StrArray *filenametab;
  PairSeqpos *filelengthtab;
  uint32_t prefixlength;
  const Seqpos *suftab;
  Encodedsequence *encseq;
  Specialcharinfo specialcharinfo;
  Alphabet *alpha;
} Suffixarray;

#endif
