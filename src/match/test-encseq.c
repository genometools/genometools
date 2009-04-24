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

#include "core/error.h"
#include "core/array.h"
#include "core/sequence_buffer_fasta.h"
#include "core/sequence_buffer_plain.h"
#include "core/progressbar.h"
#include "spacedef.h"
#include "readmode-def.h"
#include "encseq-def.h"
#include "stamp.h"

#include "arrcmp.pr"

static void runscanatpostrial(const Encodedsequence *encseq,
                              Encodedsequencescanstate *esr,
                              Readmode readmode,Seqpos startpos)
{
  Seqpos pos, totallength;
  GtUchar ccra, ccsr;

  totallength = getencseqtotallength(encseq);
  initEncodedsequencescanstate(esr,encseq,readmode,startpos);
  for (pos=startpos; pos < totallength; pos++)
  {
    ccra = getencodedchar(encseq,pos,readmode); /* Random access */
    ccsr = sequentialgetencodedchar(encseq,esr,pos,readmode);
    if (ccra != ccsr)
    {
      fprintf(stderr,"startpos = " FormatSeqpos
                     " access=%s, mode=%s: position=" FormatSeqpos
                     ": random access (correct) = %u != %u = "
                     " sequential read (wrong)\n",
                     PRINTSeqposcast(startpos),
                     encseqaccessname(encseq),
                     showreadmode(readmode),
                     PRINTSeqposcast(pos),
                     (unsigned int) ccra,
                     (unsigned int) ccsr);
      exit(EXIT_FAILURE); /* programming error */
    }
  }
}

static void testscanatpos(const Encodedsequence *encseq,
                          Readmode readmode,
                          unsigned long scantrials)
{
  Encodedsequencescanstate *esr = NULL;
  Seqpos startpos, totallength;
  unsigned long trial;

  totallength = getencseqtotallength(encseq);
  srand48(42349421);
  esr = newEncodedsequencescanstate();
  runscanatpostrial(encseq,esr,readmode,0);
  runscanatpostrial(encseq,esr,readmode,totallength-1);
  for (trial = 0; trial < scantrials; trial++)
  {
    startpos = (Seqpos) (drand48() * (double) totallength);
    printf("trial %lu at " FormatSeqpos "\n",trial,PRINTSeqposcast(startpos));
    runscanatpostrial(encseq,esr,readmode,startpos);
  }
  freeEncodedsequencescanstate(&esr);
}

static void testmulticharactercompare(const Encodedsequence *encseq,
                                      Readmode readmode,
                                      unsigned long multicharcmptrials)
{
  Encodedsequencescanstate *esr1, *esr2;
  Seqpos pos1, pos2, totallength;
  unsigned long trial;
  bool fwd = ISDIRREVERSE(readmode) ? false : true,
       complement = ISDIRCOMPLEMENT(readmode) ? true : false;

  esr1 = newEncodedsequencescanstate();
  esr2 = newEncodedsequencescanstate();
  totallength = getencseqtotallength(encseq);
  srand48(42349421);
  (void) multicharactercompare_withtest(encseq,fwd,complement,esr1,0,esr2,0);
  (void) multicharactercompare_withtest(encseq,fwd,complement,esr1,0,esr2,
                                        totallength-1);
  (void) multicharactercompare_withtest(encseq,fwd,complement,esr1,
                                        totallength-1,esr2,0);
  (void) multicharactercompare_withtest(encseq,fwd,complement,esr1,
                                        totallength-1,esr2,totallength-1);
  for (trial = 0; trial < multicharcmptrials; trial++)
  {
    pos1 = (Seqpos) (drand48() * (double) totallength);
    pos2 = (Seqpos) (drand48() * (double) totallength);
    (void) multicharactercompare_withtest(encseq,fwd,complement,
                                          esr1,pos1,esr2,pos2);
  }
  freeEncodedsequencescanstate(&esr1);
  freeEncodedsequencescanstate(&esr2);
}

static int testfullscan(const GtStrArray *filenametab,
                        const Encodedsequence *encseq,
                        Readmode readmode,
                        GtError *err)
{
  Seqpos pos, totallength;
  GtUchar ccscan = 0, ccra, ccsr;
  GtSequenceBuffer *fb = NULL;
  int retval;
  bool haserr = false;
  Encodedsequencescanstate *esr;
  unsigned long long fullscanpbar = 0;

  gt_error_check(err);
  totallength = getencseqtotallength(encseq);
  gt_progressbar_start(&fullscanpbar,(unsigned long long) totallength);
  if (filenametab != NULL)
  {
    fb = gt_sequence_buffer_new_guess_type((GtStrArray*) filenametab, err);
    if (!fb)
      haserr = true;
    if (!haserr)
      gt_sequence_buffer_set_symbolmap(fb, getencseqAlphabetsymbolmap(encseq));
  }
  if (!haserr) {
    esr = newEncodedsequencescanstate();
    initEncodedsequencescanstate(esr,encseq,readmode,0);
    for (pos=0; /* Nothing */; pos++)
    {
      if (filenametab != NULL && readmode == Forwardmode)
      {
        retval = gt_sequence_buffer_next(fb,&ccscan,err);
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
          gt_error_set(err,"access=%s, position=" FormatSeqpos
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
      ccsr = sequentialgetencodedchar(encseq,esr,pos,readmode);
      if (ccra != ccsr)
      {
        gt_error_set(err,"access=%s, mode=%s: position=" FormatSeqpos
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
    gt_progressbar_stop();
  }
  if (!haserr)
  {
    if (pos != totallength)
    {
      gt_error_set(err,"sequence length must be " FormatSeqpos " but is "
                         FormatSeqpos,totallength,pos);
      haserr = true;
    }
  }
  freeEncodedsequencescanstate(&esr);
  gt_sequence_buffer_delete(fb);
  return haserr ? -1 : 0;
}

int testencodedsequence(const GtStrArray *filenametab,
                        const Encodedsequence *encseq,
                        Readmode readmode,
                        unsigned long scantrials,
                        unsigned long multicharcmptrials,
                        GtError *err)
{
  bool fwd = ISDIRREVERSE(readmode) ? false : true,
       complement = ISDIRCOMPLEMENT(readmode) ? true : false;

  if (hasfastspecialrangeenumerator(encseq))
  {
    checkextractunitatpos(encseq,fwd,complement);
    if (multicharcmptrials > 0)
    {
      testmulticharactercompare(encseq,readmode,multicharcmptrials);
    }
  }
  if (!complement)
  {
    checkextractspecialbits(encseq,fwd);
  }
  if (scantrials > 0)
  {
    testscanatpos(encseq,readmode,scantrials);
  }
  return testfullscan(filenametab,encseq,readmode,err);
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
  GtArray *rangesforward, *rangesbackward;
  bool haserr = false;
  Specialrangeiterator *sri;
  Sequencerange range;

  if (!hasspecialranges(encseq))
  {
    return 0;
  }
  rangesforward = gt_array_new(sizeof (Sequencerange));
  rangesbackward = gt_array_new(sizeof (Sequencerange));

  sri = newspecialrangeiterator(encseq,true);
  while (nextspecialrangeiterator(&range,sri))
  {
    gt_array_add(rangesforward,range);
  }
  freespecialrangeiterator(&sri);
  sri = newspecialrangeiterator(encseq,false);
  while (nextspecialrangeiterator(&range,sri))
  {
    gt_array_add(rangesbackward,range);
  }
  freespecialrangeiterator(&sri);
  gt_array_reverse(rangesbackward);
  if (!haserr)
  {
    if (array_compare(rangesforward,rangesbackward,
                      compareSequencerange) != 0)
    {
      exit(EXIT_FAILURE); /* programming error */
    }
  }
  gt_array_delete(rangesforward);
  gt_array_delete(rangesbackward);
  return haserr ? - 1 : 0;
}
