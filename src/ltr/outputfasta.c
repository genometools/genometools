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
#include "match/echoseq.pr"

#include "ltrharvest-opt.h"
#include "repeattypes.h"

#include "match/pos2seqnum.pr"

/*
 The datetype Fastaoutinfo aggregates the values needed to show a fasta file.
 */

#define FASTASEPARATOR '>'

typedef struct
{
  Sequentialsuffixarrayreader *ssar;
  const Encodedsequence *encseq; /* encoded sequence */
  const Alphabet *alpha;         /* the alphabet */
  const Uchar *characters;       /* for visible characters */
  const char *destab;            /* pointer on descriptions */
  unsigned long *descendtab;     /* positions of desc-separators */
  Seqpos totallength;            /* totallength of encseq */
  unsigned long numofdbsequences; /* num of sequences in suffix array */
  unsigned int linewidth;        /* the line width to show the alignment */
  Seqpos *markpos;     /* positions of SEPARATOR symbols */
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
                         const char *desc,
                         const Alphabet *alpha,
                         const Encodedsequence *encseq,
                         Readmode readmode,
                         Seqpos start,
                         Seqpos wlen,
                         unsigned long width)
{
  Seqpos offset;

  assert(width > 0);
  if (desc == NULL)
  {
    fprintf(fpout,">");
  } else
  {
    fprintf(fpout,">%s",desc);
  }

  fprintf(info->formatout," (dbseq-nr");
  if (info->showseqnum)
  {
    fprintf(info->formatout," %lu) ", seqnum);
  }
  if (seqnum == 0)
  {
    offset = 0;
  }
  else
  {
    offset = info->markpos[seqnum-1]+1;
  }
  fprintf(info->formatout,"[" FormatSeqpos "," FormatSeqpos "]\n",
                       /* increase by one for output */
                       PRINTSeqposcast(start - offset + 1),
                       /* increase by one for output */
                       PRINTSeqposcast(start - offset + wlen));
  encseq2symbolstring(fpout,
                      alpha,
                      encseq,
                      readmode,
                      start,
                      wlen,
                      width);
}

static int showpredictionfastasequence(Fastaoutinfo *info, Seqpos startpos,
                                       Seqpos len,
                                       GT_UNUSED GtStr *str_indexfilename,
                                       GtError *err)
{
  unsigned long i, desclen;
  const char *desptr;
  char *destab_seqnum;
  unsigned long seqnum =
                  getrecordnumSeqpos(info->markpos, info->numofdbsequences,
                                     info->totallength, startpos, err);

  /* if there are sequence descriptions */
  desptr = retriesequencedescription(&desclen,
                                     info->destab,
                                     info->descendtab,
                                     seqnum);

  ALLOCASSIGNSPACE(destab_seqnum, NULL, char, desclen);
  for (i=0; i < desclen; i++)
  {
    destab_seqnum[i] = desptr[i];
  }

  myencseq2symbolstring(info, info->formatout,
                      seqnum, destab_seqnum,
                      info->alpha, info->encseq,
                      Forwardmode, startpos,
                      len,
                      60UL);
  FREESPACE(destab_seqnum);
  return 0;
}

/*
 The following function processes all predicted LTR elements with
 the function from the apply pointer.
 */
static int overallpredictionsequences(const LTRharvestoptions *lo,
    bool innerregion,
    void *applyinfo,
    int(*apply)(Fastaoutinfo *,Seqpos, Seqpos, GtStr*, GtError *err),
    GtError *err)
{
  unsigned long i;
  Seqpos start,
         end;
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
      }
      else
      {
        start = boundaries->leftLTR_5;
        end = boundaries->rightLTR_3;
      }
      if (apply(applyinfo,
               innerregion ? boundaries->leftLTR_3 + 1: boundaries->leftLTR_5,
               end - start + 1, lo->str_indexname, err) != 0)
      {
        return -1;
      }
    }
  }
  return 0;
}

/*
 The following function prints sequences in multiple FASTA format to the
 specified output.
*/
int showpredictionsmultiplefasta(const LTRharvestoptions *lo,
                                 Seqpos *markpos,
                                 bool innerregion,
                                 unsigned int linewidth,
                                 Sequentialsuffixarrayreader *ssar,
                                 bool showseqnum,
                                 GtError *err)
{
  Fastaoutinfo fastaoutinfo;
  FILE *formatout = NULL;
  unsigned long *descendtab = NULL,
                destablength;
  const char *destab = NULL;
  int had_err;

  formatout = gt_fa_xfopen(innerregion
                           ? gt_str_get(lo->str_fastaoutputfilenameinnerregion)
                           : gt_str_get(lo->str_fastaoutputfilename),
                           "w");

  fastaoutinfo.ssar = ssar;
  fastaoutinfo.encseq = encseqSequentialsuffixarrayreader(ssar);
  fastaoutinfo.alpha = alphabetSequentialsuffixarrayreader(ssar);
  fastaoutinfo.characters = getcharactersAlphabet(fastaoutinfo.alpha);
  fastaoutinfo.totallength = getencseqtotallength(
                               encseqSequentialsuffixarrayreader(ssar));
  fastaoutinfo.numofdbsequences =
                      numofdbsequencesSequentialsuffixarrayreader(ssar);
  fastaoutinfo.linewidth = linewidth;
  fastaoutinfo.showseqnum = showseqnum;
  fastaoutinfo.formatout = formatout;
  fastaoutinfo.markpos = markpos;

  destablength = destablengthSequentialsuffixarrayreader(ssar);
  destab = destabSequentialsuffixarrayreader(ssar);
  descendtab = calcdescendpositions(destab,
                                    destablength,
                                    fastaoutinfo.numofdbsequences);
  fastaoutinfo.destab = destab;
  fastaoutinfo.descendtab = descendtab;

  had_err = overallpredictionsequences(lo, innerregion, &fastaoutinfo,
                                       showpredictionfastasequence, err);

  gt_free(descendtab);
  gt_fa_xfclose(formatout);

  return 0;
}
