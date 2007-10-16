/*
  Copyright (c) 2007 David Ellinghaus <dellinghaus@zbh.uni-hamburg.de>
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

#include "libgtcore/env.h"
#include "libgtcore/str.h"
#include "libgtmatch/alphadef.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/spacedef.h"
#include "libgtmatch/encseq-def.h"
#include "libgtmatch/esa-seqread.h"
#include "libgtmatch/esa-seqread.pr"
#include "libgtmatch/pos2seqnum.pr"
#include "libgtmatch/echoseq.pr"

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
                         unsigned long wlen,
                         unsigned long width,
			 Env *env)
{
  unsigned long j;
  Seqpos idx, lastpos;
  Uchar currentchar;
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
                       PRINTSeqposcast(start - offset + (Seqpos)wlen) );

  lastpos = start + wlen - 1;
  for (idx = start, j = 0; ; idx++)
  {
    currentchar = getencodedchar(encseq,idx,readmode);
    if (currentchar == (Uchar) SEPARATOR)
    {
      fprintf(fpout,"\n>\n");
      j = 0;
    } else
    {
      echoprettysymbol(fpout,alpha,currentchar);
    }
    if (idx == lastpos)
    {
      fprintf(fpout,"\n");
      break;
    }
    if (currentchar != (Uchar) SEPARATOR)
    {
      j++;
      if (j >= width)
      {
        fprintf(fpout,"\n");
        j = 0;
      }
    }
  }
}

static int showpredictionfastasequence(Fastaoutinfo *info, Seqpos startpos,
                    Seqpos len, /*@unused@*/Str *str_indexfilename, Env *env)
{
  Seqpos i;
  unsigned long desclen;
  const char *desptr;
  char *destab_seqnum;
  unsigned long seqnum =
                  getrecordnumSeqpos(info->markpos, info->numofdbsequences,
		                     info->totallength, startpos, env);

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
		      (unsigned long)len,
		      60,
		      env);
  FREESPACE(destab_seqnum);

  return 0;
}

/*
 The following function processes all predicted LTR elements with
 the function from the apply pointer.
 */
static int overallpredictionsequences(const LTRharvestoptions *lo,
    Sequentialsuffixarrayreader *ssar,
    bool innerregion,
    void *applyinfo,
    int(*apply)(Fastaoutinfo *,Seqpos, Seqpos, Str*, Env *env),
    Env *env)
{
  unsigned int i;
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
	       end - start + 1, lo->str_indexname, env) != 0)
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
		                 Env *env)
{
  Fastaoutinfo fastaoutinfo;
  FILE *formatout = NULL;
  unsigned long *descendtab = NULL,
                destablength;
  const char *destab = NULL;
  int had_err;

  formatout = env_fa_xfopen(env,
                            innerregion
                            ? str_get(lo->str_fastaoutputfilenameinnerregion)
                            : str_get(lo->str_fastaoutputfilename),
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
                                    fastaoutinfo.numofdbsequences,
                                    env);
  fastaoutinfo.destab = destab;
  fastaoutinfo.descendtab = descendtab;

  had_err = overallpredictionsequences(lo, ssar, innerregion, &fastaoutinfo,
				       showpredictionfastasequence, env);

  env_ma_free(descendtab, env);
  env_fa_xfclose(formatout, env);

  return 0;
}
