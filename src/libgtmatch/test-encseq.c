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
#include "libgtcore/fastabuffer.h"
#include "spacedef.h"
#include "readmode-def.h"
#include "encseq-def.h"
#include "stamp.h"

#include "arrcmp.pr"

static int testscanatpos(const Encodedsequence *encseq,
                         Readmode readmode,
                         unsigned long trials,
                         Env *env)
{
  Encodedsequencescanstate *esr = NULL;
  Seqpos pos, startpos, totallength;
  unsigned long trial;
  Uchar ccra, ccsr;
  bool haserr = false;

  env_error_check(env);
  totallength = getencseqtotallength(encseq);
  srand48(42349421);
  for (trial = 0; !haserr && trial < trials; trial++)
  {
    startpos = (Seqpos) (drand48() * (double) totallength);
    printf("trial %lu at " FormatSeqpos "\n",trial,PRINTSeqposcast(startpos));
    esr = initEncodedsequencescanstate(encseq,readmode,startpos,env);
    for (pos=startpos; !haserr && pos < totallength; pos++)
    {
      ccra = getencodedchar(encseq,pos,readmode);
      ccsr = sequentialgetencodedchar(encseq,esr,pos);
      if (ccra != ccsr)
      {
        env_error_set(env,"startpos = " FormatSeqpos
                          " access=%s, mode=%s: position=" FormatSeqpos
                          ": random access (getencodedchar) = %u != "
                          " %u = sequential read (sequentialgetencodedchar)",
                          startpos,
                          encseqaccessname(encseq),
                          showreadmode(readmode),
                          pos,
                          (unsigned int) ccra,
                          (unsigned int) ccsr);
        haserr = true;
      }
    }
    freeEncodedsequencescanstate(&esr,env);
  }
  return haserr ? -1 : 0;
}

static int testfullscan(const StrArray *filenametab,
                        const Encodedsequence *encseq,
                        Readmode readmode,
                        const Uchar *symbolmap,
                        Env *env)
{
  Seqpos pos, totallength;
  Uchar ccscan = 0, ccra, ccsr;
  FastaBuffer *fb = NULL;
  int retval;
  bool haserr = false;
  Encodedsequencescanstate *esr;

  env_error_check(env);
  totallength = getencseqtotallength(encseq);
  if (filenametab != NULL)
  {
    fb = fastabuffer_new(filenametab,
                         symbolmap,
                         false,
                         NULL,
                         NULL,
                         NULL,
                         env);
  }
  esr = initEncodedsequencescanstate(encseq,readmode,0,env);
  for (pos=0; /* Nothing */; pos++)
  {
    if (filenametab != NULL && readmode == Forwardmode)
    {
      retval = fastabuffer_next(fb,&ccscan,env);
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
  fastabuffer_delete(fb, env);
  return haserr ? -1 : 0;
}

int testencodedsequence(const StrArray *filenametab,
                        const Encodedsequence *encseq,
                        Readmode readmode,
                        const Uchar *symbolmap,
                        unsigned long trials,
                        Env *env)
{
  bool haserr = false;

  if (trials > 0)
  {
    if (testscanatpos(encseq,
                     Forwardmode /* XXX readmode */,
                     trials,
                     env) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    if (testfullscan(filenametab,encseq,readmode,symbolmap,env) != 0)
    {
      haserr = true;
    }
  }
  return haserr ? -1 : 0;
}

static void makeerrormsg(const Sequencerange *vala,const Sequencerange *valb,
                         const char *cmpflag)
{
  fprintf(stderr,
                "(" FormatSeqpos "," FormatSeqpos
                ") %s (" FormatSeqpos "," FormatSeqpos
                ")\n",
                PRINTSeqposcast(vala->leftpos),
                PRINTSeqposcast(vala->rightpos),
                cmpflag,
                PRINTSeqposcast(valb->leftpos),
                PRINTSeqposcast(valb->rightpos));
}

static int compareSequencerange(const void *a,const void *b)
{
  const Sequencerange *vala, *valb;

  vala = (Sequencerange *) a;
  valb = (Sequencerange *) b;
  if (vala->leftpos < valb->leftpos)
  {
    makeerrormsg(vala,valb,"<");
    return -1;
  }
  if (vala->leftpos > valb->leftpos)
  {
    makeerrormsg(vala,valb,">");
    return 1;
  }
  if (vala->rightpos < valb->rightpos)
  {
    makeerrormsg(vala,valb,"<");
    return -1;
  }
  if (vala->rightpos > valb->rightpos)
  {
    makeerrormsg(vala,valb,">");
    return 1;
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
    if (array_compare(rangesforward,rangesbackward,
                      compareSequencerange) != 0)
    {
      exit(EXIT_FAILURE);
    }
  }
  array_delete(rangesforward,env);
  array_delete(rangesbackward,env);
  return haserr ? - 1 : 0;
}
