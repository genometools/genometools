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
#include "core/unused_api.h"
#include "match/alphadef.h"
#include "match/sarr-def.h"
#include "match/spacedef.h"
#include "match/encseq-def.h"
#include "match/esa-seqread.h"
#include "match/echoseq.h"
#include "ltrharvest-opt.h"
#include "repeattypes.h"

/*
 The datetype Fastaoutinfo aggregates the values needed to show a fasta file.
 */

#define FASTASEPARATOR '>'

typedef struct
{
  Sequentialsuffixarrayreader *ssar;
  const Encodedsequence *encseq; /* encoded sequence */
  unsigned long *descendtab;     /* positions of desc-separators */
  unsigned long linewidth;        /* the line width to show the alignment */
  bool showseqnum;     /* with or without the sequence number */
  FILE *formatout;     /* file pointer to show the alignment */
} Fastaoutinfo;

/*
 The following function processes one predicted LTR element, i.e.
 a FASTA entry is created from the predicted LTR element.
 */

static void myencseq2symbolstring(Fastaoutinfo *info,
                                  FILE *fpout,
                                  unsigned long seqnum,
                                  Seqpos offset,
                                  const char *desc,
                                  unsigned long desclength,
                                  const Encodedsequence *encseq,
                                  Readmode readmode,
                                  Seqpos start,
                                  Seqpos wlen,
                                  unsigned long linewidth)
{
  gt_assert(linewidth > 0);
  if (desc == NULL)
  {
    fprintf(fpout,">");
  } else
  {
    fprintf(fpout,">%*.*s",(int) desclength,(int) desclength,desc);
  }

  fprintf(info->formatout," (dbseq-nr");
  if (info->showseqnum)
  {
    fprintf(info->formatout," %lu) ", seqnum);
  }
  fprintf(info->formatout,"[" FormatSeqpos "," FormatSeqpos "]\n",
                       /* increase by one for output */
                       PRINTSeqposcast(start - offset + 1),
                       /* increase by one for output */
                       PRINTSeqposcast(start - offset + wlen));
  encseq2symbolstring(fpout,
                      encseq,
                      readmode,
                      start,
                      wlen,
                      linewidth);
}

static void showpredictionfastasequence(Fastaoutinfo *info, Seqpos startpos,
                                        Seqpos len)
{
  unsigned long desclen;
  const char *desptr;
  Seqinfo seqinfo;
  unsigned long seqnum = getencseqfrompos2seqnum(info->encseq,startpos);

  /* if there are sequence descriptions */
  desptr = retrievesequencedescription(&desclen,
                                       info->encseq,
                                       info->descendtab,
                                       seqnum);

  getencseqSeqinfo(&seqinfo,info->encseq,seqnum);
  myencseq2symbolstring(info,
                        info->formatout,
                        seqnum, 
                        seqinfo.seqstartpos,
                        desptr,
                        desclen,
                        info->encseq,
                        Forwardmode, 
                        startpos,
                        len,
                        info->linewidth);
}

/*
 The following function processes all predicted LTR elements with
 the function from the apply pointer.
 */

static void overallpredictionsequences(const LTRharvestoptions *lo,
                                       bool innerregion,
                                       void *applyinfo,
                                       void(*apply)(Fastaoutinfo *,Seqpos,
                                                    Seqpos))
{
  unsigned long i;
  Seqpos start, end;
  LTRboundaries *boundaries;

  for (i = 0; i < lo->arrayLTRboundaries.nextfreeLTRboundaries; i++)
  {
    boundaries = &(lo->arrayLTRboundaries.spaceLTRboundaries[i]);
    if ( !boundaries->skipped )
    {
      if (innerregion)
      {
        start = boundaries->leftLTR_3 + 1;
        end = boundaries->rightLTR_5 - 1;
      } else
      {
        start = boundaries->leftLTR_5;
        end = boundaries->rightLTR_3;
      }
      apply(applyinfo,
            innerregion ? boundaries->leftLTR_3 + 1: boundaries->leftLTR_5,
            end - start + 1);
    }
  }
}

/*
 The following function prints sequences in multiple FASTA format to the
 specified output.
*/

int showpredictionsmultiplefasta(const LTRharvestoptions *lo,
                                 bool innerregion,
                                 unsigned int linewidth,
                                 Sequentialsuffixarrayreader *ssar,
                                 bool showseqnum,
                                 GtError *err)
{
  Fastaoutinfo fastaoutinfo;
  FILE *formatout = NULL;
  bool had_err = false;

  formatout = gt_fa_fopen(innerregion
                          ? gt_str_get(lo->str_fastaoutputfilenameinnerregion)
                          : gt_str_get(lo->str_fastaoutputfilename),
                          "w",err);
  if (formatout == NULL)
  {
    had_err = true;
  }

  if (!had_err)
  {
    fastaoutinfo.ssar = ssar;
    fastaoutinfo.encseq = encseqSequentialsuffixarrayreader(ssar);
    fastaoutinfo.linewidth = (unsigned long) linewidth;
    fastaoutinfo.showseqnum = showseqnum;
    fastaoutinfo.formatout = formatout;
    fastaoutinfo.descendtab = calcdescendpositions(fastaoutinfo.encseq);
    overallpredictionsequences(lo, innerregion, &fastaoutinfo,
                               showpredictionfastasequence);
    gt_free(fastaoutinfo.descendtab);
  }
  gt_fa_xfclose(formatout);
  return had_err ? -1 : 0;
}
