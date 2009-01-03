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

#include "core/fa.h"
#include "core/str.h"
#include "match/esa-seqread.h"
#include "match/echoseq.h"
#include "ltr/ltrharvest-opt.h"
#include "ltr/repeattypes.h"

void printgff3format(LTRharvestoptions *lo, Sequentialsuffixarrayreader *ssar,
                     const Seqpos *markpos)
{
  LTRboundaries *boundaries;
  Seqpos contiglen,
         offset;
  unsigned long seqnum,
                i,
                contignumber,
       idcounterRepregion = 0,
       idcounterRetrotrans = 0,
       idcounterLTR = 0,
       idcounterTSD = 0,
       idcounterMotif = 0;

  unsigned long *descendtab = NULL,
                desclen;
  const char *desptr = NULL;

  unsigned long numofdbsequences;
  const Encodedsequence *encseq;
  Seqpos totallength;
  FILE *fp;

  encseq = encseqSequentialsuffixarrayreader(ssar);
  numofdbsequences = getencseqnumofdbsequences(encseq);
  totallength = getencseqtotallength(encseq);

  fp = gt_fa_xfopen(gt_str_get(lo->str_gff3filename), "w");

  /* for getting descriptions */
  descendtab = calcdescendpositions(encseq);

  if (lo->arrayLTRboundaries.nextfreeLTRboundaries == 0)
  {
    /* no LTR-pairs predicted */
  }
  else
  {
    fprintf(fp, "##gff-version 3\n");
    /* print output sorted by contignumber */
    for (seqnum = 0; seqnum < numofdbsequences; seqnum++)
    {
      /* contig is first sequence, and only one sequence in multiseq */
      if ( seqnum == 0 && numofdbsequences == 1UL)
      {
        contiglen = totallength;
      }
      else
      {
        /* first sequence and more than one sequence in suffixarray */
        if ( seqnum == 0)
        {
          contiglen = markpos[seqnum];
        }
        else
        {
          /* last sequence in suffixarray */
          if (seqnum == numofdbsequences - 1)
          {
            contiglen = totallength - 1 - markpos[seqnum-1];
          }
          else
          {
            contiglen = markpos[seqnum] - markpos[seqnum-1] - 1;
          }
        }
      }
      fprintf(fp, "##sequence-region seq%lu 1 " FormatSeqpos "\n",
                  seqnum, PRINTSeqposcast(contiglen));
      /* write description of sequence */
      fprintf(fp, "# ");
      desptr = retrievesequencedescription(&desclen,
                                           encseq,
                                           descendtab,
                                           seqnum);
      for (i=0; i < desclen; i++)
      {
        fprintf(fp, "%c", desptr[i]);
      }
      fprintf(fp, "\n");

      for (i = 0; i < lo->arrayLTRboundaries.nextfreeLTRboundaries; i++)
      {
        boundaries = &(lo->arrayLTRboundaries.spaceLTRboundaries[i]);
        contignumber = boundaries->contignumber;
        if ( (!boundaries->skipped) && contignumber == seqnum)
        {
          if ( contignumber == 0)
          {
            offset = 0;
          }
          else
          {
            offset = markpos[contignumber-1]+1;
          }

          /* repeat-region */
          fprintf(fp, "seq%lu\tLTRharvest\trepeat_region\t" FormatSeqpos
                      "\t" FormatSeqpos "\t" ".\t?\t.\tID=RepeatReg%lu\n",
              contignumber,
              /* increase boundary position by one for output */
              PRINTSeqposcast(boundaries->leftLTR_5 -offset + 1
                          - boundaries->lenleftTSD),
              PRINTSeqposcast(boundaries->rightLTR_3 -offset + 1
                          + boundaries->lenrightTSD),
              idcounterRepregion++ );

          /* LTR retrotransposon */
          fprintf(fp, "seq%lu\tLTRharvest\tLTR_retrotransposon\t"
                      FormatSeqpos "\t" FormatSeqpos "\t"
                      ".\t?\t.\tID=LTRret%lu;Parent=RepeatReg%lu\n",
              contignumber,
              /* increase boundary position by one for output */
              PRINTSeqposcast(boundaries->leftLTR_5 -offset + 1),
              PRINTSeqposcast(boundaries->rightLTR_3 -offset  + 1),
              idcounterRetrotrans++,
              idcounterRepregion-1 );

          /* LTRs */
          fprintf(fp, "seq%lu\tLTRharvest\tlong_terminal_repeat\t"
                      FormatSeqpos "\t" FormatSeqpos "\t"
                      ".\t?\t.\tID=LTR%lu;Parent=LTRret%lu\n",
              contignumber,
              /* increase boundary position by one for output */
              PRINTSeqposcast(boundaries->leftLTR_5 -offset + 1),
              PRINTSeqposcast(boundaries->leftLTR_3 -offset + 1),
              idcounterLTR++,
              idcounterRetrotrans-1 );
          fprintf(fp, "seq%lu\tLTRharvest\tlong_terminal_repeat\t"
                      FormatSeqpos "\t" FormatSeqpos "\t"
                      ".\t?\t.\tID=LTR%lu;Parent=LTRret%lu\n",
              contignumber,
              /* increase boundary position by one for output */
              PRINTSeqposcast(boundaries->rightLTR_5 -offset + 1),
              PRINTSeqposcast(boundaries->rightLTR_3 -offset + 1),
              idcounterLTR++,
              idcounterRetrotrans-1 );

          if (lo->minlengthTSD > 1U)
          {
            /* TSDs */
            fprintf(fp, "seq%lu\tLTRharvest\ttarget_site_duplication\t"
                        FormatSeqpos "\t" FormatSeqpos "\t"
                        ".\t?\t.\tID=TSD%lu;Parent=RepeatReg%lu\n",
                contignumber,
                /* increase boundary position by one for output */
                PRINTSeqposcast(boundaries->leftLTR_5 -offset + 1
                          - boundaries->lenleftTSD),
                PRINTSeqposcast(boundaries->leftLTR_5 -offset),
                idcounterTSD++,
                idcounterRepregion-1 );

            fprintf(fp, "seq%lu\tLTRharvest\ttarget_site_duplication\t"
                        FormatSeqpos "\t" FormatSeqpos "\t"
                        ".\t?\t.\tID=TSD%lu;Parent=RepeatReg%lu\n",
                contignumber,
                /* increase boundary position by one for output */
                PRINTSeqposcast(boundaries->rightLTR_3 -offset + 2),
                PRINTSeqposcast(boundaries->rightLTR_3 -offset + 1
                          + boundaries->lenrightTSD),
                idcounterTSD++,
                idcounterRepregion-1 );
          }

          if (lo->motif.allowedmismatches < 4U)
          {
            fprintf(fp, "seq%lu\tLTRharvest\tinverted_repeat\t"
                        FormatSeqpos "\t" FormatSeqpos "\t"
                        ".\t?\t.\tID=Motif%lu;Parent=RepeatReg%lu\n",
                contignumber,
                /* increase boundary position by one for output */
                PRINTSeqposcast(boundaries->leftLTR_5 -offset + 1),
                PRINTSeqposcast(boundaries->leftLTR_5 -offset + 2),
                idcounterMotif++,
                idcounterRepregion-1 );
            fprintf(fp, "seq%lu\tLTRharvest\tinverted_repeat\t"
                        FormatSeqpos "\t" FormatSeqpos "\t"
                        ".\t?\t.\tID=Motif%lu;Parent=RepeatReg%lu\n",
                contignumber,
                /* increase boundary position by one for output */
                PRINTSeqposcast(boundaries->leftLTR_3 -offset),
                PRINTSeqposcast(boundaries->leftLTR_3 -offset + 1),
                idcounterMotif++,
                idcounterRepregion-1 );

            fprintf(fp, "seq%lu\tLTRharvest\tinverted_repeat\t"
                        FormatSeqpos "\t" FormatSeqpos "\t"
                        ".\t?\t.\tID=Motif%lu;Parent=RepeatReg%lu\n",
                contignumber,
                /* increase boundary position by one for output */
                PRINTSeqposcast(boundaries->rightLTR_5 -offset + 1),
                PRINTSeqposcast(boundaries->rightLTR_5 -offset + 2),
                idcounterMotif++,
                idcounterRepregion-1 );
            fprintf(fp, "seq%lu\tLTRharvest\tinverted_repeat\t"
                        FormatSeqpos "\t" FormatSeqpos "\t"
                        ".\t?\t.\tID=Motif%lu;Parent=RepeatReg%lu\n",
                contignumber,
                /* increase boundary position by one for output */
                PRINTSeqposcast(boundaries->rightLTR_3 -offset),
                PRINTSeqposcast(boundaries->rightLTR_3 -offset + 1),
                idcounterMotif++,
                idcounterRepregion-1 );
          }
        }
      }
    }
  }
  gt_free(descendtab);
  gt_fa_xfclose(fp);
}
