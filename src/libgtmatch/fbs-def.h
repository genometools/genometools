/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef FBS_DEF_H
#define FBS_DEF_H
#include <stdio.h>
#include <stdbool.h>
#include <inttypes.h>
#include "libgtcore/strarray.h"
#include "symboldef.h"
#include "filelength-def.h"
#include "arraydef.h"
#include "seqdesc.h"

#define FILEBUFFERSIZE 65536

typedef struct
{
  unsigned int filenum;
  uint64_t linenum;
  unsigned long nextread,
                nextfree;
  bool indesc,
       firstoverallseq,
       firstseqinfile,
       complete,
       nextfile;
  Sequencedescription *sequencedescription;
  GenFile *inputstream;
  Uchar bufspace[FILEBUFFERSIZE];
  uint64_t lastspeciallength;
  Filelengthvalues *filelengthtab;
  const StrArray *filenametab;
  const Uchar *symbolmap;
  bool plainformat;
  unsigned long *characterdistribution;
} Fastabufferstate;

Fastabufferstate*  initformatbufferstate(const StrArray *filenametab,
                                         const Uchar *symbolmap,
                                         bool plainformat,
                                         Filelengthvalues **filelengthtab,
                                         Sequencedescription
                                         *sequencedescription,
                                         unsigned long *characterdistribution,
                                         Env *env);

int advanceformatbufferstate(Fastabufferstate *fbs,Env *env);

static inline int readnextUchar(Uchar *val,Fastabufferstate *fbs,Env *env)
{
  if (fbs->nextread >= fbs->nextfree)
  {
    if (fbs->complete)
    {
      return 0;
    }
    if (advanceformatbufferstate(fbs,env) != 0)
    {
      return -1;
    }
    fbs->nextread = 0;
    if (fbs->nextfree == 0)
    {
      return 0;
    }
  }
  *val = fbs->bufspace[fbs->nextread++];
  return 1;
}

void fastabufferstate_delete(Fastabufferstate*, Env*);

#endif
