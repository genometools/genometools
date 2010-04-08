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
#include "core/encodedsequence.h"
#include "ltrharvest-opt.h"
#include "outputstd.h"

/*
   The following function prints the predicted LTR retrotransposon
   results to stdout.
 */

/* CAUTION: For output the positions will be incremened by one,
 *          since normally in annotation first base in sequence
 *          is position 1 instead of 0.
 */

static void printlongheader(const LTRharvestoptions *lo)
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
}

static void producelongutput(const LTRharvestoptions *lo,
                             const LTRboundaries *boundaries,
                             const GtEncodedsequence *encseq,
                             unsigned long offset)
{
  const GtUchar *characters = gt_alphabet_characters(
                                           gt_encodedsequence_alphabet(encseq));

  printf("%lu  ",
      boundaries->leftLTR_5 -offset + 1);
  printf("%lu  ",
      boundaries->rightLTR_3 -offset  + 1);
  printf("%lu  ",
      boundaries->rightLTR_3 - boundaries->leftLTR_5 + 1);
  printf("%lu  ",
      boundaries->leftLTR_5 -offset  + 1);
  printf("%lu  ",
      boundaries->leftLTR_3 -offset  + 1);
  printf("%lu  ",
      boundaries->leftLTR_3 - boundaries->leftLTR_5 + 1);
  if (lo->minlengthTSD > 1U)
  {
    unsigned long j;

    for (j = 0; j < boundaries->lenleftTSD; j++)
    {
      printf("%c",(char) characters[gt_encodedsequence_get_encoded_char(encseq,
                                                   boundaries->leftLTR_5 -
                                                   boundaries->lenleftTSD + j,
                                                   GT_READMODE_FORWARD)]);
    }
    printf("  " "%lu  ",
           boundaries->lenleftTSD);
  }
  if (lo->motif.allowedmismatches < 4U)
  {
    printf("%c%c..%c%c  ",
        /* Random access */
        (char) characters[gt_encodedsequence_get_encoded_char(encseq,
                       boundaries->leftLTR_5,
                       GT_READMODE_FORWARD)],
        /* Random access */
        (char) characters[gt_encodedsequence_get_encoded_char(encseq,
                       boundaries->leftLTR_5+1,
                       GT_READMODE_FORWARD)],
        /* Random access */
        (char) characters[gt_encodedsequence_get_encoded_char(encseq,
                       boundaries->leftLTR_3-1,
                       GT_READMODE_FORWARD)],
        /* Random access */
        (char) characters[gt_encodedsequence_get_encoded_char(encseq,
                       boundaries->leftLTR_3,
                       GT_READMODE_FORWARD)] );
  }
  /* increase by 1 */
  printf("%lu  ",
      boundaries->rightLTR_5 -offset + 1);
  /* increase by 1 */
  printf("%lu  ",boundaries->rightLTR_3 -offset + 1);
  printf("%lu  ",boundaries->rightLTR_3
                                           - boundaries->rightLTR_5 + 1);
  if (lo->minlengthTSD > 1U)
  {
    unsigned long j;

    for (j = 0; j < boundaries->lenrightTSD; j++)
    {
      printf("%c", (char) characters[gt_encodedsequence_get_encoded_char(encseq,
                                                    boundaries->rightLTR_3+j+1,
                                                    GT_READMODE_FORWARD)]);
    }
    printf("  " "%lu  ",boundaries->lenrightTSD);
  }
  if (lo->motif.allowedmismatches < 4U)
  {
    printf("%c%c..%c%c",
        /* Random access */
        (char) characters[gt_encodedsequence_get_encoded_char(encseq,
                       boundaries->rightLTR_5,
                       GT_READMODE_FORWARD)],
        /* Random access */
        (char) characters[gt_encodedsequence_get_encoded_char(encseq,
                       boundaries->rightLTR_5+1,
                       GT_READMODE_FORWARD)],
        /* Random access */
        (char) characters[gt_encodedsequence_get_encoded_char(encseq,
                       boundaries->rightLTR_3-1,
                       GT_READMODE_FORWARD)],
        /* Random access */
        (char) characters[gt_encodedsequence_get_encoded_char(encseq,
                       boundaries->rightLTR_3,
                       GT_READMODE_FORWARD)] );
  }
  /* print similarity */
  printf("  %.2f", boundaries->similarity);
  /* print sequence number */
  printf("  %lu\n", boundaries->contignumber);
}

static void printshortheader(void)
{
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
}

static void produceshortoutput(const LTRboundaries *boundaries,
                               unsigned long offset)
{

  /* increase positions by 1 */
  printf("%lu  ",
      boundaries->leftLTR_5 -offset + 1);
  printf("%lu  ",
       boundaries->rightLTR_3 -offset + 1);
  printf("%lu  ",
      boundaries->rightLTR_3
          - boundaries->leftLTR_5 + 1);
  printf("%lu  ",
      boundaries->leftLTR_5 -offset + 1);
  printf("%lu  ",
      boundaries->leftLTR_3 -offset + 1);
  printf("%lu  ",
      boundaries->leftLTR_3
          - boundaries->leftLTR_5 + 1);
  printf("%lu  ",
      boundaries->rightLTR_5 -offset + 1);
  printf("%lu  ",
      boundaries->rightLTR_3 -offset + 1);
  printf("%lu  ",
      boundaries->rightLTR_3
          - boundaries->rightLTR_5 + 1);
  /* print similarity */
  printf("%.2f  ", boundaries->similarity);
  /* print sequence number */
  printf("%lu\n", boundaries->contignumber);
}

void gt_showinfoiffoundfullLTRs(const LTRharvestoptions *lo,
                             const LTRboundaries **bdptrtab,
                             unsigned long numofboundaries,
                             const GtEncodedsequence *encseq)
{
  GtSeqinfo seqinfo;
  unsigned long i;

  if (lo->longoutput)
  {
    if (numofboundaries == 0)
    {
      printf("No full LTR-pair predicted.\n");
    } else
    {
      printlongheader(lo);
      /* print output sorted by contignumber*/
      for (i = 0; i<numofboundaries; i++)
      {
        gt_encodedsequence_seqinfo(encseq,&seqinfo,bdptrtab[i]->contignumber);
        producelongutput(lo,
                         bdptrtab[i],
                         encseq,
                         seqinfo.seqstartpos);
      }
    }
  } else
  {
    if (numofboundaries > 0)
    {
      printshortheader();
      /* print output sorted by contignumber*/
      for (i = 0; i<numofboundaries; i++)
      {
        gt_encodedsequence_seqinfo(encseq,&seqinfo,bdptrtab[i]->contignumber);
        produceshortoutput(bdptrtab[i],seqinfo.seqstartpos);
      }
    }
  }
}
