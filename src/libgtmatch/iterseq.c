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

#include <inttypes.h>
#include "libgtcore/strarray.h"
#include "symboldef.h"
#include "arraydef.h"
#include "chardef.h"
#include "fbs-def.h"
#include "seqdesc.h"
#include "gqueue-def.h"
#include "format64.h"
#include "iterseq.h"

#include "fbsadv.pr"
#include "genericqueue.pr"

#include "readnextUchar.gen"

int overallquerysequences(int(*processsequence)(void *,
                                                uint64_t,
                                                const Uchar *,
                                                unsigned long,
                                                const char *,
                                                Env *),
                          void *info,
                          ArrayUchar *sequencebuffer,
                          const StrArray *filenametab,
                          Sequencedescription *sequencedescription,
                          const Uchar *symbolmap,
                          Env *env)
{
  Fastabufferstate fbs;
  Uchar charcode;
  int retval;
  uint64_t unitnum = 0;
  bool haserr = false;
  char *desc;

  env_error_check(env);
  initformatbufferstate(&fbs,
                        filenametab,
                        symbolmap,
                        false,
                        NULL,
                        sequencedescription,
                        env);
  sequencebuffer->nextfreeUchar = 0;
  while (true)
  {
    retval = readnextUchar(&charcode,&fbs,env);
    if (retval < 0)
    {
      haserr = true;
      break;
    }
    if (retval == 0)
    {
      break;
    }
    if (charcode == (Uchar) SEPARATOR)
    {
      if (sequencebuffer->nextfreeUchar == 0)
      {
        env_error_set(env,"sequence " Formatuint64_t " is empty",
                      PRINTuint64_tcast(unitnum));
        haserr = true;
        break;
      }
      desc = dequeuegeneric(sequencedescription->descptr,env);
      if (desc == NULL)
      {
        haserr = true;
        break;
      }
      if (processsequence(info,
                          unitnum,
                          sequencebuffer->spaceUchar,
                          sequencebuffer->nextfreeUchar,
                          desc,
                          env) != 0)
      {
        haserr = true;
        FREESPACE(desc);
        break;
      }
      FREESPACE(desc);
      sequencebuffer->nextfreeUchar = 0;
      unitnum++;
    } else
    {
      STOREINARRAY(sequencebuffer,Uchar,1024,charcode);
    }
  }
  if (!haserr && sequencebuffer->nextfreeUchar > 0)
  {
    desc = dequeuegeneric(sequencedescription->descptr,env);
    if (desc == NULL)
    {
      haserr = true;
    }
    if (!haserr)
    {
      if (processsequence(info,
                          unitnum,
                          sequencebuffer->spaceUchar,
                          sequencebuffer->nextfreeUchar,
                          desc,
                          env) != 0)
      {
        haserr = true;
      }
    }
    FREESPACE(desc);
    sequencebuffer->nextfreeUchar = 0;
  }
  return haserr ? -1 : 0;
}

 struct Scansequenceiterator
{
  Fastabufferstate fbs;
  const StrArray *filenametab;
  const Uchar *symbolmap;
  Sequencedescription sequencedescription;
  ArrayUchar sequencebuffer;
  uint64_t unitnum;
  bool exhausted;
};

Scansequenceiterator *newScansequenceiterator(const StrArray *filenametab,
                                              const Uchar *symbolmap,
                                              Env *env)
{
  Scansequenceiterator *sseqit;

  ALLOCASSIGNSPACE(sseqit,NULL,Scansequenceiterator,1);
  INITARRAY(&sseqit->sequencebuffer,Uchar);
  INITARRAY(&sseqit->sequencedescription.headerbuffer,char);
  sseqit->sequencedescription.descptr = emptyqueuegeneric(env);
  initformatbufferstate(&sseqit->fbs,
                        filenametab,
                        symbolmap,
                        false,
                        NULL,
                        &sseqit->sequencedescription,
                        env);
  sseqit->sequencebuffer.nextfreeUchar = 0;
  sseqit->exhausted = false;
  sseqit->unitnum = 0;
  return sseqit;
}

void freeScansequenceiterator(Scansequenceiterator **sseqit,Env *env)
{
  wrapqueuegeneric(true,&(*sseqit)->sequencedescription.descptr,env);
  FREEARRAY(&(*sseqit)->sequencebuffer,Uchar);
  FREEARRAY(&(*sseqit)->sequencedescription.headerbuffer,char);
  FREESPACE(*sseqit);
}

int nextScansequenceiterator(const Uchar **sequence,
                             unsigned long *len,
                             char **desc,
                             Scansequenceiterator *sseqit,
                             Env *env)
{
  Uchar charcode;
  int retval;
  bool haserr = false, foundseq = false;

  if(sseqit->exhausted)
  {
    return 0;
  }
  while (true)
  {
    retval = readnextUchar(&charcode,&sseqit->fbs,env);
    if (retval < 0)
    {
      haserr = true;
      break;
    }
    if (retval == 0)
    {
      sseqit->exhausted = true;
      break;
    }
    if (charcode == (Uchar) SEPARATOR)
    {
      if (sseqit->sequencebuffer.nextfreeUchar == 0)
      {
        env_error_set(env,"sequence " Formatuint64_t " is empty",
                      PRINTuint64_tcast(sseqit->unitnum));
        haserr = true;
        break;
      }
      *desc = dequeuegeneric(sseqit->sequencedescription.descptr,env);
      if (*desc == NULL)
      {
        haserr = true;
        break;
      }
      *sequence = sseqit->sequencebuffer.spaceUchar;
      *len = sseqit->sequencebuffer.nextfreeUchar;
      foundseq = true;
      sseqit->sequencebuffer.nextfreeUchar = 0;
      foundseq = true;
      sseqit->unitnum++;
      break;
    } else
    {
      STOREINARRAY(&sseqit->sequencebuffer,Uchar,1024,charcode);
    }
  }
  if (!haserr && sseqit->sequencebuffer.nextfreeUchar > 0)
  {
    *desc = dequeuegeneric(sseqit->sequencedescription.descptr,env);
    if (*desc == NULL)
    {
      haserr = true;
    }
    if (!haserr)
    {
      *sequence = sseqit->sequencebuffer.spaceUchar;
      *len = sseqit->sequencebuffer.nextfreeUchar;
      foundseq = true;
    }
    sseqit->sequencebuffer.nextfreeUchar = 0;
  }
  if(haserr)
  {
    return -1;
  }
  if(foundseq)
  {
    return 1;
  }
  return 0;
}
