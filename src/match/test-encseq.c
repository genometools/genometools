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
                              GtReadmode readmode,unsigned long startpos)
{
  unsigned long pos, totallength;
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
      fprintf(stderr,"startpos = %lu"
                     " access=%s, mode=%s: position=%lu"
                     ": random access (correct) = %u != %u = "
                     " sequential read (wrong)\n",
                     startpos,
                     encseqaccessname(encseq),
                     gt_readmode_show(readmode),
                     pos,
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
  unsigned long startpos, totallength;
  unsigned long trial;

  totallength = gt_encodedsequence_total_length(encseq);
  srand48(42349421);
  esr = gt_encodedsequence_scanstate_new();
  runscanatpostrial(encseq,esr,readmode,0);
  runscanatpostrial(encseq,esr,readmode,totallength-1);
  for (trial = 0; trial < scantrials; trial++)
  {
    startpos = (unsigned long) (drand48() * (double) totallength);
    printf("trial %lu at %lu\n",trial,startpos);
    runscanatpostrial(encseq,esr,readmode,startpos);
  }
  gt_encodedsequence_scanstate_delete(esr);
}

static void testmulticharactercompare(const GtEncodedsequence *encseq,
                                      GtReadmode readmode,
                                      unsigned long multicharcmptrials)
{
  GtEncodedsequenceScanstate *esr1, *esr2;
  unsigned long pos1, pos2, totallength;
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
    pos1 = (unsigned long) (drand48() * (double) totallength);
    pos2 = (unsigned long) (drand48() * (double) totallength);
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
  unsigned long pos, totallength;
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
          gt_error_set(err,"access=%s, position=%lu"
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
        gt_error_set(err,"access=%s, mode=%s: position=%lu"
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
      gt_error_set(err,"sequence length must be %lu but is %lu",
                       totallength,pos);
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

  if (gt_encodedsequence_has_fast_specialrangeenumerator(encseq))
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

static void makeerrormsg(const GtRange *vala,
                         const GtRange *valb,
                         const char *cmpflag)
{
  fprintf(stderr,
                "(%lu,%lu) %s (%lu,%lu)\n",
                vala->start,
                vala->end,
                cmpflag,
                valb->start,
                valb->end);
}

static int compareGtRange(const void *a,const void *b)
{
  const GtRange *vala, *valb;

  vala = (GtRange *) a;
  valb = (GtRange *) b;
  if (vala->start < valb->start)
  {
    makeerrormsg(vala,valb,"<");
    return -1;
  }
  if (vala->start > valb->start)
  {
    makeerrormsg(vala,valb,">");
    return 1;
  }
  if (vala->end < valb->end)
  {
    makeerrormsg(vala,valb,"<");
    return -1;
  }
  if (vala->end > valb->end)
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
  GtSpecialrangeiterator *sri;
  GtRange range;

  if (!gt_encodedsequence_has_specialranges(encseq))
  {
    return 0;
  }
  rangesforward = gt_array_new(sizeof (GtRange));
  rangesbackward = gt_array_new(sizeof (GtRange));

  sri = gt_specialrangeiterator_new(encseq,true);
  while (gt_specialrangeiterator_next(sri,&range))
  {
    gt_array_add(rangesforward,range);
  }
  gt_specialrangeiterator_delete(sri);
  sri = gt_specialrangeiterator_new(encseq,false);
  while (gt_specialrangeiterator_next(sri,&range))
  {
    gt_array_add(rangesbackward,range);
  }
  gt_specialrangeiterator_delete(sri);
  gt_array_reverse(rangesbackward);
  if (!haserr)
  {
    if (array_compare(rangesforward,rangesbackward,
                      compareGtRange) != 0)
    {
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
  gt_array_delete(rangesforward);
  gt_array_delete(rangesbackward);
  return haserr ? - 1 : 0;
}
