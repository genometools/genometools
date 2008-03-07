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

#include "libgtcore/error.h"
#include "libgtcore/array.h"
#include "libgtcore/fastabuffer.h"
#include "libgtcore/progressbar.h"
#include "spacedef.h"
#include "readmode-def.h"
#include "encseq-def.h"
#include "stamp.h"

#include "arrcmp.pr"

static void testscanatpos(const Encodedsequence *encseq,
                          Readmode readmode,
                          unsigned long trials)
{
  Encodedsequencescanstate *esr = NULL;
  Seqpos pos, startpos, totallength;
  unsigned long trial;
  Uchar ccra, ccsr;

  totallength = getencseqtotallength(encseq);
  srand48(42349421);
  for (trial = 0; trial < trials; trial++)
  {
    startpos = (Seqpos) (drand48() * (double) totallength);
    printf("trial %lu at " FormatSeqpos "\n",trial,PRINTSeqposcast(startpos));
    esr = newEncodedsequencescanstate();
    initEncodedsequencescanstate(esr,encseq,readmode,startpos);
    for (pos=startpos; pos < totallength; pos++)
    {
      ccra = getencodedchar(encseq,pos,readmode); /* Random access */
      ccsr = sequentialgetencodedchar(encseq,esr,pos);
      if (ccra != ccsr)
      {
        fprintf(stderr,"startpos = " FormatSeqpos
                       " access=%s, mode=%s: position=" FormatSeqpos
                       ": random access (correct) = %u != %u = "
                       " sequential read (wrong)",
                       PRINTSeqposcast(startpos),
                       encseqaccessname(encseq),
                       showreadmode(readmode),
                       PRINTSeqposcast(pos),
                       (unsigned int) ccra,
                       (unsigned int) ccsr);
        exit(EXIT_FAILURE); /* programming error */
      }
    }
    freeEncodedsequencescanstate(&esr);
  }
}

static void testmulticharactercompare(const Encodedsequence *encseq,
                                      unsigned long trials)
{
  Encodedsequencescanstate *esr1, *esr2;
  Seqpos pos1, pos2, totallength;
  unsigned long trial;
  int ret1, ret2;

  esr1 = newEncodedsequencescanstate();
  esr2 = newEncodedsequencescanstate();
  totallength = getencseqtotallength(encseq);
  srand48(42349421);
  for (trial = 0; trial < trials; trial++)
  {
    pos1 = (Seqpos) (drand48() * (double) totallength);
    pos2 = (Seqpos) (drand48() * (double) totallength);
    ret1 = multicharactercompare(encseq,esr1,pos1,esr2,pos2);
    ret2 = multicharactercompare_bruteforce(encseq,pos1,pos2);
    if (ret1 != ret2)
    {
      fprintf(stderr,"pos1=" FormatSeqpos ", pos2=" FormatSeqpos "\n",
              PRINTSeqposcast(pos1),PRINTSeqposcast(pos2));
      fprintf(stderr,"ret1=%d, ret2=%d\n",ret1,ret2);
      exit(EXIT_FAILURE); /* programming error */
    }
  }
  freeEncodedsequencescanstate(&esr1);
  freeEncodedsequencescanstate(&esr2);
}

static int testfullscan(const StrArray *filenametab,
                        const Encodedsequence *encseq,
                        Readmode readmode,
                        const Uchar *symbolmap,
                        Error *err)
{
  Seqpos pos, totallength;
  Uchar ccscan = 0, ccra, ccsr;
  FastaBuffer *fb = NULL;
  int retval;
  bool haserr = false;
  Encodedsequencescanstate *esr;
  unsigned long long fullscanpbar = 0;

  error_check(err);
  totallength = getencseqtotallength(encseq);
  progressbar_start(&fullscanpbar,(unsigned long long) totallength);
  if (filenametab != NULL)
  {
    fb = fastabuffer_new(filenametab,
                         symbolmap,
                         false,
                         NULL,
                         NULL,
                         NULL);
  }
  esr = newEncodedsequencescanstate();
  initEncodedsequencescanstate(esr,encseq,readmode,0);
  for (pos=0; /* Nothing */; pos++)
  {
    if (filenametab != NULL && readmode == Forwardmode)
    {
      retval = fastabuffer_next(fb,&ccscan,err);
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
    ccra = getencodedchar(encseq,pos,readmode); /* Random access */
    if (filenametab != NULL && readmode == Forwardmode)
    {
      if (ccscan != ccra)
      {
        error_set(err,"access=%s, position=" FormatSeqpos
                          ": scan (readnextchar) = %u != "
                          "%u = random access",
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
      error_set(err,"access=%s, mode=%s: position=" FormatSeqpos
                        ": random access = %u != %u = sequential read",
                        encseqaccessname(encseq),
                        showreadmode(readmode),
                        pos,
                        (unsigned int) ccra,
                        (unsigned int) ccsr);
      haserr = true;
      break;
    }
    fullscanpbar++;
  }
  progressbar_stop();
  if (!haserr)
  {
    if (pos != totallength)
    {
      error_set(err,"sequence length must be " FormatSeqpos " but is "
                         FormatSeqpos,totallength,pos);
      haserr = true;
    }
  }
  freeEncodedsequencescanstate(&esr);
  fastabuffer_delete(fb);
  return haserr ? -1 : 0;
}

int testencodedsequence(const StrArray *filenametab,
                        const Encodedsequence *encseq,
                        Readmode readmode,
                        const Uchar *symbolmap,
                        unsigned long trials,
                        Error *err)
{
  if (hasfastspecialrangeenumerator(encseq))
  {
    checkextractunitatpos(encseq);
    if (trials > 0)
    {
      testmulticharactercompare(encseq,trials);
    }
  }
  if (trials > 0)
  {
    testscanatpos(encseq,readmode,trials);
  }
  return testfullscan(filenametab,encseq,readmode,symbolmap,err);
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

int checkspecialrangesfast(const Encodedsequence *encseq)
{
  Array *rangesforward, *rangesbackward;
  bool haserr = false;
  Specialrangeiterator *sri;
  Sequencerange range;

  if (!hasspecialranges(encseq))
  {
    return 0;
  }
  rangesforward = array_new(sizeof (Sequencerange));
  rangesbackward = array_new(sizeof (Sequencerange));

  sri = newspecialrangeiterator(encseq,true);
  while (nextspecialrangeiterator(&range,sri))
  {
    array_add(rangesforward,range);
  }
  freespecialrangeiterator(&sri);
  sri = newspecialrangeiterator(encseq,false);
  while (nextspecialrangeiterator(&range,sri))
  {
    array_add(rangesbackward,range);
  }
  freespecialrangeiterator(&sri);
  array_reverse(rangesbackward);
  if (!haserr)
  {
    if (array_compare(rangesforward,rangesbackward,
                      compareSequencerange) != 0)
    {
      exit(EXIT_FAILURE);
    }
  }
  array_delete(rangesforward);
  array_delete(rangesbackward);
  return haserr ? - 1 : 0;
}
