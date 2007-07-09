/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FBS_DEF_H
#define FBS_DEF_H
#include <stdio.h>
#include <stdbool.h>
#include <zlib.h>
#include "libgtcore/strarray.h"
#include "types.h"

#define FILEBUFFERSIZE 65536

typedef struct
{
  bool isgzippedstream;
  union
  {
    FILE *fopenstream;
    gzFile gzippedstream;
  } stream;
} Genericstream;

typedef struct
{
  uint32_t filenum,
           linenum,
           nextread,
           nextfree;
  bool indesc,
       firstoverallseq,
       firstseqinfile,
       complete,
       nextfile;
  Genericstream inputstream;
  Uchar bufspace[FILEBUFFERSIZE];
  Seqpos totaloffset;
  Seqpos lastspeciallength;
  PairSeqpos *filelengthtab;
  const StrArray *filenametab;
  const Uchar *symbolmap;
} Fastabufferstate;

#endif
