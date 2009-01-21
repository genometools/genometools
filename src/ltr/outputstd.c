/*
  Copyright (c) 2007 David Ellinghaus <d.ellinghaus@ikmb.uni-kiel.de>
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

#include "core/symboldef.h"
#include "match/sarr-def.h"
#include "match/encseq-def.h"
#include "match/readmode-def.h"
#include "match/intcode-def.h"
#include "match/spacedef.h"
#include "ltrharvest-opt.h"

/*
   The following function prints the predicted LTR retrotransposon
   results to stdout.
 */

/* CAUTION: For output the positions will be incremened by one,
 *          since normally in annotation first base in sequence
 *          is position 1 instead of 0.
 */

static void producelongutput(const LTRharvestoptions *lo,
                             const LTRboundaries *boundaries,
                             const Encodedsequence *encseq,
                             Seqpos offset)
{
  const Uchar *characters = getencseqAlphabetcharacters(encseq);

  printf(FormatSeqpos "  ",
      PRINTSeqposcast(boundaries->leftLTR_5 -offset + (Seqpos)1));
  printf(FormatSeqpos "  ",
      PRINTSeqposcast(boundaries->rightLTR_3 -offset  + (Seqpos)1));
  printf(FormatSeqpos "  ",
      PRINTSeqposcast((boundaries->rightLTR_3 - boundaries->leftLTR_5
          + 1)));
  printf(FormatSeqpos "  ",
      PRINTSeqposcast(boundaries->leftLTR_5 -offset  + (Seqpos)1));
  printf(FormatSeqpos "  ",
      PRINTSeqposcast(boundaries->leftLTR_3 -offset  + (Seqpos)1));
  printf(FormatSeqpos "  ",
      PRINTSeqposcast((boundaries->leftLTR_3 - boundaries->leftLTR_5
          + (Seqpos)1)));
  if (lo->minlengthTSD > 1U)
  {
    Seqpos j;

    for (j = 0; j < boundaries->lenleftTSD; j++)
    {
      printf("%c",
       (char) characters[getencodedchar(encseq,/* XXX */
          boundaries->leftLTR_5 - boundaries->lenleftTSD + j,
          Forwardmode)]);
    }
    printf("  " FormatSeqpos "  ",
           PRINTSeqposcast(boundaries->lenleftTSD));
  }
  if (lo->motif.allowedmismatches < 4U)
  {
    printf("%c%c..%c%c  ",
        (char) characters[getencodedchar(encseq,/* Random access */
                       boundaries->leftLTR_5,
                       Forwardmode)],
        (char) characters[getencodedchar(encseq,/* Random access */
                       boundaries->leftLTR_5+(Seqpos)1,
                       Forwardmode)],
        (char) characters[getencodedchar(encseq,/* Random access */
                       boundaries->leftLTR_3-(Seqpos)1,
                       Forwardmode)],
        (char) characters[getencodedchar(encseq,/* Random access */
                       boundaries->leftLTR_3,
                       Forwardmode)] );
  }
  /* increase by 1 */
  printf(FormatSeqpos "  ",
      PRINTSeqposcast(boundaries->rightLTR_5 -offset + (Seqpos)1));
  /* increase by 1 */
  printf(FormatSeqpos "  ",
      PRINTSeqposcast(boundaries->rightLTR_3 -offset + (Seqpos)1));
  printf(FormatSeqpos "  ",
      PRINTSeqposcast((boundaries->rightLTR_3
          - boundaries->rightLTR_5 + 1)));
  if (lo->minlengthTSD > 1U)
  {
    Seqpos j;

    for (j = 0; j < boundaries->lenrightTSD; j++)
    {
      printf("%c", (char) characters[getencodedchar(encseq,/* XXX */
          boundaries->rightLTR_3 + j + 1,
          Forwardmode)]);
    }
    printf("  " FormatSeqpos "  ",
           PRINTSeqposcast(boundaries->lenrightTSD));
  }
  if (lo->motif.allowedmismatches < 4U)
  {
    printf("%c%c..%c%c",
        (char) characters[getencodedchar(encseq,/* Randomaccess */
                       boundaries->rightLTR_5,
                       Forwardmode)],
        (char) characters[getencodedchar(encseq,/* Randomaccess */
                       boundaries->rightLTR_5+(Seqpos)1,
                       Forwardmode)],
        (char) characters[getencodedchar(encseq,/* Randomaccess */
                       boundaries->rightLTR_3-(Seqpos)1,
                       Forwardmode)],
        (char) characters[getencodedchar(encseq,/* Random access */
                       boundaries->rightLTR_3,/* Randomaccess */
                       Forwardmode)] );
  }
  /* print similarity */
  printf("  %.2f", boundaries->similarity);
  /* print sequence number */
  printf("  %lu\n", boundaries->contignumber);
}

static void produceshortoutput(const LTRboundaries *boundaries,Seqpos offset)
{

  /* increase positions by 1 */
  printf(FormatSeqpos "  ",
      PRINTSeqposcast( boundaries->leftLTR_5 -offset + (Seqpos)1));
  printf(FormatSeqpos "  ",
      PRINTSeqposcast( boundaries->rightLTR_3 -offset + (Seqpos)1));
  printf(FormatSeqpos "  ",
      PRINTSeqposcast( (boundaries->rightLTR_3
          - boundaries->leftLTR_5 + (Seqpos)1)));
  printf(FormatSeqpos "  ",
      PRINTSeqposcast( boundaries->leftLTR_5 -offset + (Seqpos)1));
  printf(FormatSeqpos "  ",
      PRINTSeqposcast( boundaries->leftLTR_3 -offset + (Seqpos)1));
  printf(FormatSeqpos "  ",
      PRINTSeqposcast( (boundaries->leftLTR_3
          - boundaries->leftLTR_5 + (Seqpos)1)));
  printf(FormatSeqpos "  ",
      PRINTSeqposcast( boundaries->rightLTR_5 -offset + (Seqpos)1));
  printf(FormatSeqpos "  ",
      PRINTSeqposcast( boundaries->rightLTR_3 -offset + (Seqpos)1));
  printf(FormatSeqpos "  ",
      PRINTSeqposcast( (boundaries->rightLTR_3
          - boundaries->rightLTR_5 + (Seqpos)1)));
  /* print similarity */
  printf("%.2f  ", boundaries->similarity);
  /* print sequence number */
  printf("%lu\n", boundaries->contignumber);
}

int showinfoiffoundfullLTRs(const LTRharvestoptions *lo,
                            const Sequentialsuffixarrayreader *ssar)
{
  LTRboundaries *boundaries;
  Seqinfo seqinfo;
  unsigned long numofdbsequences,
                seqnum,
                i;
  const Encodedsequence *encseq;

  encseq = encseqSequentialsuffixarrayreader(ssar);
  /* in order to get to visible dna characters */

  numofdbsequences = getencseqnumofdbsequences(encseq);

  if (lo->longoutput)
  {
    if (lo->arrayLTRboundaries.nextfreeLTRboundaries == 0)
    {
      printf("No full LTR-pair predicted.\n");
    } else
    {
      printf("# predictions are reported in the following way\n");
      printf("# s(ret) e(ret) l(ret) ");
      printf("s(lLTR) e(lLTR) l(lLTR)");
      if (lo->minlengthTSD > 1U)
      {
        printf(" TSD l(TSD)");
      }
      if (lo->motif.allowedmismatches < 4U)
      {
        printf(" m(lLTR)");
      }
      printf(" s(rLTR) e(rLTR) l(rLTR)");
      if (lo->minlengthTSD > 1U)
      {
        printf(" TSD l(TSD)");
      }
      if (lo->motif.allowedmismatches < 4U)
      {
        printf(" m(rLTR)");
      }
      printf(" sim(LTRs)");
      printf(" seq-nr");
      printf("\n# where:\n");
      printf("# s = starting position\n");
      printf("# e = ending position\n");
      printf("# l = length\n");
      if (lo->motif.allowedmismatches < 4U)
      {
        printf("# m = motif\n");
      }
      printf("# ret = LTR-retrotransposon\n");
      printf("# lLTR = left LTR\n");
      printf("# rLTR = right LTR\n");
      if (lo->minlengthTSD > 1U)
      {
        printf("# TSD = target site duplication\n");
      }
      printf("# sim = similarity\n");
      printf("# seq-nr = sequence number\n");

      /* print output sorted by contignumber*/
      for (seqnum = 0; seqnum < numofdbsequences; seqnum++)
      {
        for (i = 0; i < lo->arrayLTRboundaries.nextfreeLTRboundaries; i++)
        {
          boundaries = &(lo->arrayLTRboundaries.spaceLTRboundaries[i]);
          if ( (!boundaries->skipped) && boundaries->contignumber == seqnum)
          {
            getencseqSeqinfo(&seqinfo,encseq,boundaries->contignumber);
            producelongutput(lo,
                             boundaries,
                             encseq,
                             seqinfo.seqstartpos);
          }
        }
      }
    }
  } else
  {
    if (lo->arrayLTRboundaries.nextfreeLTRboundaries > 0)
    {
      /* print short output of full length LTR-retrotransposon(s) */
      printf("# predictions are reported in the following way\n");
      printf("# s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR)"
          " s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr \n");
      printf("# where:\n");
      printf("# s = starting position\n");
      printf("# e = ending position\n");
      printf("# l = length\n");
      printf("# ret = LTR-retrotransposon\n");
      printf("# lLTR = left LTR\n");
      printf("# rLTR = right LTR\n");
      printf("# sim = similarity\n");
      printf("# seq-nr = sequence number\n");

      /* print output sorted by contignumber*/
      for (seqnum = 0; seqnum < numofdbsequences; seqnum++)
      {
        for (i = 0; i < lo->arrayLTRboundaries.nextfreeLTRboundaries; i++)
        {
          boundaries = &(lo->arrayLTRboundaries.spaceLTRboundaries[i]);
          if ( (!boundaries->skipped) && boundaries->contignumber == seqnum)
          {
            getencseqSeqinfo(&seqinfo,encseq,boundaries->contignumber);
            produceshortoutput(boundaries,seqinfo.seqstartpos);
          }
        }
      }
    }
  }
  return 0;
}
