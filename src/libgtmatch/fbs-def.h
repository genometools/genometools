/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FBS_DEF_H
#define FBS_DEF_H
#include <stdbool.h>
#include "types.h"
#include "genstream.h"

#define FILEBUFFERSIZE 65536

typedef struct
{
  unsigned int filenum;
  unsigned int numoffiles;
  bool indesc,
       firstoverallseq,
       firstseqinfile,
       complete;
  Uint linenum;
  Genericstream inputstream;
  bool nextfile;
  Uint nextread,
       nextfree;
  Uchar bufspace[FILEBUFFERSIZE];
  Uint64 totaloffset;
  const char **filenametab;
  const Uchar *symbolmap;
  Uint lastspeciallength;
  PairUint *filelengthtab;
} Fastabufferstate;

#endif
