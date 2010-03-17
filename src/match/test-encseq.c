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
#include "core/readmode.h"
#include "core/encodedsequence.h"
#include "stamp.h"

#include "arrcmp.pr"

static void runscanatpostrial(const GtEncodedsequence *encseq,
                              GtEncodedsequenceScanstate *esr,
                              GtReadmode readmode,Seqpos startpos)
{
  Seqpos pos, totallength;
  GtUchar ccra, ccsr;

  totallength = gt_encodedsequence_total_length(encseq);
  gt_encodedsequence_scanstate_init(esr,encseq,readmode,startpos);
  for (pos=startpos; pos < totallength; pos++)
  {
    /* Random access */
    ccra = gt_encodedsequence_getencodedchar(encseq,pos,readmode);
    ccsr = gt_encodedsequence_sequentialgetencodedchar(encseq,esr,pos,readmode);
    if (ccra != ccsr)
    {
      fprintf(stderr,"startpos = " FormatSeqpos
                     " access=%s, mode=%s: position=" FormatSeqpos
                     ": random access (correct) = %u != %u = "
                     " sequential read (wrong)\n",
                     PRINTSeqposcast(startpos),
                     encseqaccessname(encseq),
                     gt_readmode_show(readmode),
                     PRINTSeqposcast(pos),
                     (unsigned int) ccra,
                     (unsigned int) ccsr);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
}

static void testscanatpos(const GtEncodedsequence *encseq,
                          GtReadmode readmode,
                          unsigned long scantrials)
{
  GtEncodedsequenceScanstate *esr = NULL;
  Seqpos startpos, totallength;
  unsigned long trial;

  totallength = gt_encodedsequence_total_length(encseq);
  srand48(42349421);
  esr = gt_encodedsequence_scanstate_new();
  runscanatpostrial(encseq,esr,readmode,0);
  runscanatpostrial(encseq,esr,readmode,totallength-1);
  for (trial = 0; trial < scantrials; trial++)
  {
    startpos = (Seqpos) (drand48() * (double) totallength);
    printf("trial %lu at " FormatSeqpos "\n",trial,PRINTSeqposcast(startpos));
    runscanatpostrial(encseq,esr,readmode,startpos);
  }
  gt_encodedsequence_scanstate_delete(esr);
}

static void testmulticharactercompare(const GtEncodedsequence *encseq,
                                      GtReadmode readmode,
                                      unsigned long multicharcmptrials)
{
  GtEncodedsequenceScanstate *esr1, *esr2;
  Seqpos pos1, pos2, totallength;
  unsigned long trial;
  bool fwd = GT_ISDIRREVERSE(readmode) ? false : true,
       complement = GT_ISDIRCOMPLEMENT(readmode) ? true : false;

  esr1 = gt_encodedsequence_scanstate_new();
  esr2 = gt_encodedsequence_scanstate_new();
  totallength = gt_encodedsequence_total_length(encseq);
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
  gt_encodedsequence_scanstate_delete(esr1);
  gt_encodedsequence_scanstate_delete(esr2);
}

static int testfullscan(const GtStrArray *filenametab,
                        const GtEncodedsequence *encseq,
                        GtReadmode readmode,
                        GtError *err)
{
  Seqpos pos, totallength;
  GtUchar ccscan = 0, ccra, ccsr;
  GtSequenceBuffer *fb = NULL;
  int retval;
  bool haserr = false;
  GtEncodedsequenceScanstate *esr;
  unsigned long long fullscanpbar = 0;

  gt_error_check(err);
  totallength = gt_encodedsequence_total_length(encseq);
  gt_progressbar_start(&fullscanpbar,(unsigned long long) totallength);
  if (filenametab != NULL)
  {
    fb = gt_sequence_buffer_new_guess_type((GtStrArray*) filenametab, err);
    if (!fb)
      haserr = true;
    if (!haserr)
      gt_sequence_buffer_set_symbolmap(fb,
                                  gt_encodedsequence_alphabetsymbolmap(encseq));
  }
  if (!haserr) {
    esr = gt_encodedsequence_scanstate_new();
    gt_encodedsequence_scanstate_init(esr,encseq,readmode,0);
    for (pos=0; /* Nothing */; pos++)
    {
      if (filenametab != NULL && readmode == GT_READMODE_FORWARD)
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
      /* Random access */
      ccra = gt_encodedsequence_getencodedchar(encseq,pos,readmode);
      if (filenametab != NULL && readmode == GT_READMODE_FORWARD)
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
      ccsr = gt_encodedsequence_sequentialgetencodedchar(encseq,esr,pos,
                                                         readmode);
      if (ccra != ccsr)
      {
        gt_error_set(err,"access=%s, mode=%s: position=" FormatSeqpos
                          ": random access = %u != %u = sequential read",
                          encseqaccessname(encseq),
                          gt_readmode_show(readmode),
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
  gt_encodedsequence_scanstate_delete(esr);
  gt_sequence_buffer_delete(fb);
  return haserr ? -1 : 0;
}

int testencodedsequence(const GtStrArray *filenametab,
                        const GtEncodedsequence *encseq,
                        GtReadmode readmode,
                        unsigned long scantrials,
                        unsigned long multicharcmptrials,
                        GtError *err)
{
  bool fwd = GT_ISDIRREVERSE(readmode) ? false : true,
       complement = GT_ISDIRCOMPLEMENT(readmode) ? true : false;

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

static void makeerrormsg(const GtSequencerange *vala,
                         const GtSequencerange *valb,
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

static int compareGtSequencerange(const void *a,const void *b)
{
  const GtSequencerange *vala, *valb;

  vala = (GtSequencerange *) a;
  valb = (GtSequencerange *) b;
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

int checkspecialrangesfast(const GtEncodedsequence *encseq)
{
  GtArray *rangesforward, *rangesbackward;
  bool haserr = false;
  Specialrangeiterator *sri;
  GtSequencerange range;

  if (!hasspecialranges(encseq))
  {
    return 0;
  }
  rangesforward = gt_array_new(sizeof (GtSequencerange));
  rangesbackward = gt_array_new(sizeof (GtSequencerange));

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
                      compareGtSequencerange) != 0)
    {
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
  gt_array_delete(rangesforward);
  gt_array_delete(rangesbackward);
  return haserr ? - 1 : 0;
}
