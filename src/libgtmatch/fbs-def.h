/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FBS_DEF_H
#define FBS_DEF_H
#include <stdbool.h>
#include "libgtcore/strarray.h"
#include "types.h"
#include "genstream.h"

#define FILEBUFFERSIZE 65536

typedef struct
{
  unsigned long filenum;
  bool indesc,
       firstoverallseq,
       firstseqinfile,
       complete;
  unsigned int linenum;
  Genericstream inputstream;
  bool nextfile;
  unsigned int nextread,
               nextfree;
  Uchar bufspace[FILEBUFFERSIZE];
  Uint64 totaloffset;
  const StrArray *filenametab;
  const Uchar *symbolmap;
  Seqpos lastspeciallength;
  PairSeqpos *filelengthtab;
} Fastabufferstate;

#endif
