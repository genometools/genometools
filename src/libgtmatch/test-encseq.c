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

#include "libgtcore/env.h"
#include "libgtcore/array.h"
#include "spacedef.h"
#include "encseq-def.h"
#include "fbs-def.h"
#include "stamp.h"

#include "sfx-readmode.pr"
#include "fbsadv.pr"
#include "readnextUchar.gen"

int testencodedsequence(const StrArray *filenametab,
                        const Encodedsequence *encseq,
                        Readmode readmode,
                        const Uchar *symbolmap,
                        Env *env)
{
  Seqpos pos, totallength;
  Uchar ccscan = 0, ccra, ccsr;
  Fastabufferstate fbs;
  int retval;
  bool haserr = false;
  Encodedsequencescanstate *esr;

  env_error_check(env);
  printf("# testencodedsequence with readmode = %s\n",showreadmode(readmode));
  totallength = getencseqtotallength(encseq);
  if (filenametab != NULL)
  {
    initformatbufferstate(&fbs,
                          filenametab,
                          symbolmap,
                          false,
                          NULL,
                          NULL,
                          env);
  }
  esr = initEncodedsequencescanstate(encseq,readmode,env);
  for (pos=0; /* Nothing */; pos++)
  {
    if (filenametab != NULL && readmode == Forwardmode)
    {
      retval = readnextUchar(&ccscan,&fbs,env);
      if (retval < 0)
      {
        haserr = true;
        break;
      }
      if (retval == 0)
      {
        break;
      }
    } else
    {
      if (pos >= totallength)
      {
        break;
      }
    }
    ccra = getencodedchar(encseq,pos,readmode);
    if (filenametab != NULL && readmode == Forwardmode)
    {
      if (ccscan != ccra)
      {
        env_error_set(env,"access=%s, position=" FormatSeqpos
                          ": scan (readnextchar) = %u != "
                          "%u = random access (getencodedchar)",
                          encseqaccessname(encseq),
                          pos,
                          (unsigned int) ccscan,
                          (unsigned int) ccra);
        haserr = true;
        break;
      }
    }
    ccsr = sequentialgetencodedchar(encseq,esr,pos);
    if (ccra != ccsr)
    {
      env_error_set(env,"access=%s, mode=%s: position=" FormatSeqpos
                        ": random access (getencodedchar) = %u != "
                        " %u = sequential read (sequentialgetencodedchar)",
                        encseqaccessname(encseq),
                        showreadmode(readmode),
                        pos,
                        (unsigned int) ccra,
                        (unsigned int) ccsr);
      haserr = true;
      break;
    }
  }
  if (!haserr)
  {
    if (pos != totallength)
    {
      env_error_set(env,"sequence length must be " FormatSeqpos " but is "
                         FormatSeqpos,totallength,pos);
      haserr = true;
    }
  }
  freeEncodedsequencescanstate(&esr,env);
  return haserr ? -1 : 0;
}

static int comparearrays(const Array *a,const Array *b,Env *env)
{
  unsigned long idx;
  Sequencerange *vala, *valb;

  if (array_size(a) != array_size(b))
  {
    env_error_set(env,"array_size(a) = %lu != %lu = "
                      "array_size(b)\n",
    (unsigned long) array_size(a),
    (unsigned long) array_size(b));
    return -1;
  }
  for (idx=0; idx<(unsigned long) array_size(a); idx++)
  {
    vala = (Sequencerange *) array_get(a,idx);
    valb = (Sequencerange *) array_get(b,idx);
    if (vala->leftpos != valb->leftpos || vala->rightpos != valb->rightpos)
    {
      env_error_set(env,
                    "a[%lu] = (" FormatSeqpos "," FormatSeqpos
                      ") != (" FormatSeqpos "," FormatSeqpos
                      ") = b[%lu]",
                      idx,
                      PRINTSeqposcast(vala->leftpos),
                      PRINTSeqposcast(vala->rightpos),
                      PRINTSeqposcast(valb->leftpos),
                      PRINTSeqposcast(valb->rightpos),
                      idx);
      return -1;
    }
  }
  return 0;
}

int checkspecialrangesfast(const Encodedsequence *encseq,Env *env)
{
  Array *rangesforward, *rangesbackward;
  bool haserr = false;
  Specialrangeiterator *sri;
  Sequencerange range;

  env_error_check(env);
  if (!hasspecialranges(encseq))
  {
    return 0;
  }
  rangesforward = array_new(sizeof (Sequencerange),env);
  rangesbackward = array_new(sizeof (Sequencerange),env);

  sri = newspecialrangeiterator(encseq,true,env);
  while (nextspecialrangeiterator(&range,sri))
  {
    array_add(rangesforward,range,env);
  }
  freespecialrangeiterator(&sri,env);
  sri = newspecialrangeiterator(encseq,false,env);
  while (nextspecialrangeiterator(&range,sri))
  {
    array_add(rangesbackward,range,env);
  }
  freespecialrangeiterator(&sri,env);
  array_reverse(rangesbackward,env);
  if (!haserr)
  {
    printf("# checkspecialrangesfast(%lu ranges)\n",
             (unsigned long) array_size(rangesforward));
    if (comparearrays(rangesforward,rangesbackward,env) != 0)
    {
      haserr = true;
    }
  }
  array_delete(rangesforward,env);
  array_delete(rangesbackward,env);
  return haserr ? - 1 : 0;
}
