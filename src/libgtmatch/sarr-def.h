#ifndef SARR_DEF_H
#define SARR_DEF_H

#include "types.h"
#include "alphadef.h"
#include "encseq-def.h"

typedef struct
{
  Uint numofdbsequences;
  char **filenametab;
  PairUint *filelengthtab;
  unsigned int numoffiles;
  unsigned int prefixlength;
  const Uint *suftab;
  Encodedsequence *encseq;
  Specialcharinfo specialcharinfo;
  Alphabet *alpha;
} Suffixarray;

#endif
