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

#include "core/error.h"
#include "core/fa.h"
#include "core/str.h"
#include "core/encodedsequence.h"
#include "match/echoseq.h"
#include "outputfasta.h"
#include "ltrharvest-opt.h"

/*
 The datetype Fastaoutinfo aggregates the values needed to show a fasta file.
 */

#define FASTASEPARATOR '>'

typedef struct
{
  const GtEncodedsequence *encseq; /* encoded sequence */
  unsigned long linewidth;        /* the line width to show the alignment */
  bool showseqnum;     /* with or without the sequence number */
  FILE *formatout;     /* file pointer to show the alignment */
} Fastaoutinfo;

/*
 The following function processes one predicted LTR element, i.e.
 a FASTA entry is created from the predicted LTR element.
 */

static void myencseq2symbolstring(Fastaoutinfo *fastainfo,
                                  unsigned long seqnum,
                                  Seqpos offset,
                                  const char *desc,
                                  unsigned long desclength,
                                  GtReadmode readmode,
                                  Seqpos start,
                                  Seqpos wlen)
{
  gt_assert(fastainfo->linewidth > 0);
  if (desc == NULL)
  {
    fprintf(fastainfo->formatout,">");
  } else
  {
    fprintf(fastainfo->formatout,">%*.*s",
            (int) desclength,(int) desclength,desc);
  }

  fprintf(fastainfo->formatout," (dbseq-nr");
  if (fastainfo->showseqnum)
  {
    fprintf(fastainfo->formatout," %lu) ", seqnum);
  }
  fprintf(fastainfo->formatout,"[" FormatSeqpos "," FormatSeqpos "]\n",
                       /* increase by one for output */
                       PRINTSeqposcast(start - offset + 1),
                       /* increase by one for output */
                       PRINTSeqposcast(start - offset + wlen));
  encseq2symbolstring(fastainfo->formatout,
                      fastainfo->encseq,
                      readmode,
                      start,
                      wlen,
                      fastainfo->linewidth);
}

static void showpredictionfastasequence(Fastaoutinfo *fastainfo,
                                        Seqpos startpos,
                                        Seqpos len)
{
  unsigned long desclen;
  const char *desptr;
  GtSeqinfo seqinfo;
  unsigned long seqnum = getencseqfrompos2seqnum(fastainfo->encseq,startpos);

  desptr = gt_encodedsequence_description(fastainfo->encseq,&desclen,seqnum);
  gt_encodedsequence_seqinfo(fastainfo->encseq,&seqinfo,seqnum);
  myencseq2symbolstring(fastainfo,
                        seqnum,
                        seqinfo.seqstartpos,
                        desptr,
                        desclen,
                        GT_READMODE_FORWARD,
                        startpos,
                        len);
}

/*
 The following function processes all predicted LTR elements with
 the function from the apply pointer.
 */

static void overallpredictionsequences(const LTRboundaries **bdptrtab,
                                       unsigned long numofboundaries,
                                       bool innerregion,
                                       Fastaoutinfo *fastainfo)
{
  unsigned long i;
  Seqpos start, end;

  for (i = 0; i < numofboundaries; i++)
  {
    if (innerregion)
    {
      start = bdptrtab[i]->leftLTR_3 + 1;
      end = bdptrtab[i]->rightLTR_5 - 1;
    } else
    {
      start = bdptrtab[i]->leftLTR_5;
      end = bdptrtab[i]->rightLTR_3;
    }
    showpredictionfastasequence(fastainfo,
                                innerregion ? bdptrtab[i]->leftLTR_3 + 1
                                            : bdptrtab[i]->leftLTR_5,
                                end - start + 1);
  }
}

/*
 The following function prints sequences in multiple FASTA format to the
 specified output.
*/

int showpredictionsmultiplefasta(const LTRharvestoptions *lo,
                                 const LTRboundaries **bdptrtab,
                                 unsigned long numofboundaries,
                                 bool innerregion,
                                 unsigned int linewidth,
                                 bool showseqnum,
                                 GtError *err)
{
  FILE *formatout = NULL;
  bool had_err = false;

  formatout = gt_fa_fopen(innerregion
                          ? gt_str_get(lo->str_fastaoutputfilenameinnerregion)
                          : gt_str_get(lo->str_fastaoutputfilename),
                          "w",err);
  if (formatout == NULL)
  {
    had_err = true;
  } else
  {
    Fastaoutinfo fastaoutinfo;

    fastaoutinfo.encseq = lo->repeatinfo.encseq;
    fastaoutinfo.linewidth = (unsigned long) linewidth;
    fastaoutinfo.showseqnum = showseqnum;
    fastaoutinfo.formatout = formatout;
    overallpredictionsequences(bdptrtab, numofboundaries, innerregion,
                               &fastaoutinfo);
  }
  gt_fa_xfclose(formatout);
  return had_err ? -1 : 0;
}
