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
  Filelengthvalues *filelengthtab;
  uint32_t prefixlength;
  const Seqpos *suftab;
  const Uchar *lcptab;
  const Largelcpvalue *llvtab;
  const Uchar *bwttab;
  DefinedSeqpos numoflargelcpvalues;
  Encodedsequence *encseq;
  DefinedSeqpos longest;
  Specialcharinfo specialcharinfo;
  Alphabet *alpha;
} Suffixarray;

#endif
