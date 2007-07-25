/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FBS_DEF_H
#define FBS_DEF_H
#include <stdio.h>
#include <stdbool.h>
#include <inttypes.h>
#include <zlib.h>
#include "libgtcore/strarray.h"
#include "symboldef.h"
#include "filelength-def.h"
#include "arraydef.h"
#include "seqdesc.h"

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
  uint32_t filenum;
  uint64_t linenum;
  unsigned long nextread,
                nextfree;
  bool indesc,
       firstoverallseq,
       firstseqinfile,
       complete,
       nextfile;
  Sequencedescription *sequencedescription;
  Genericstream inputstream;
  Uchar bufspace[FILEBUFFERSIZE];
  uint64_t lastspeciallength;
  Filelengthvalues *filelengthtab;
  const StrArray *filenametab;
  const Uchar *symbolmap;
  bool plainformat;
} Fastabufferstate;

#endif
