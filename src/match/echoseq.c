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

#include "core/alphabet.h"
#include "core/assert_api.h"
#include "core/chardef.h"
#include "core/error.h"
#include "core/seqiterator.h"
#include "spacedef.h"
#include "encseq-def.h"

void symbolstring2fasta(FILE *fpout,
                        const char *desc,
                        const GtAlphabet *alpha,
                        const GtUchar *w,
                        unsigned long wlen,
                        unsigned long width)
{
  unsigned long i, j;
  GtUchar currentchar;

  gt_assert(width > 0);
  if (desc == NULL)
  {
    fprintf(fpout,">\n");
  } else
  {
    fprintf(fpout,">%s\n",desc);
  }
  for (i = 0, j = 0; ; i++)
  {
    currentchar = w[i];
    if (currentchar == (GtUchar) SEPARATOR)
    {
      fprintf(fpout,"\n>\n");
      j = 0;
    } else
    {
      gt_alphabet_echo_pretty_symbol(alpha,fpout,currentchar);
    }
    if (i == wlen - 1)
    {
      fprintf(fpout,"\n");
      break;
    }
    if (currentchar != (GtUchar) SEPARATOR)
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

void encseq2symbolstring(FILE *fpout,
                         const Encodedsequence *encseq,
                         Readmode readmode,
                         Seqpos start,
                         Seqpos wlen,
                         unsigned long width)
{
  unsigned long j;
  Seqpos idx, lastpos;
  GtUchar currentchar;
  Encodedsequencescanstate *esr;
  const GtAlphabet *alpha;

  esr = newEncodedsequencescanstate();
  initEncodedsequencescanstate(esr,encseq,readmode,start);
  gt_assert(width > 0);
  lastpos = start + wlen - 1;
  alpha = getencseqAlphabet(encseq);
  for (idx = start, j = 0; /* Nothing */ ; idx++)
  {
    currentchar = sequentialgetencodedchar(encseq,esr,idx,readmode);
    if (currentchar == (GtUchar) SEPARATOR)
    {
      fprintf(fpout,"\n>\n");
      j = 0;
    } else
    {
      gt_alphabet_echo_pretty_symbol(alpha,fpout,currentchar);
    }
    if (idx == lastpos)
    {
      fprintf(fpout,"\n");
      break;
    }
    if (currentchar != (GtUchar) SEPARATOR)
    {
     j++;
     if (j >= width)
      {
        fprintf(fpout,"\n");
        j = 0;
      }
    }
  }
  freeEncodedsequencescanstate(&esr);
}

void fprintfencseq(FILE *fpout,
                   const Encodedsequence *encseq,
                   Seqpos start,
                   Seqpos wlen)
{
  Seqpos idx;
  GtUchar currentchar;
  const GtAlphabet *alpha;

  alpha = getencseqAlphabet(encseq);
  for (idx = start; idx < start + wlen; idx++)
  {
    currentchar = getencodedchar(encseq,idx,Forwardmode);
    gt_assert(ISNOTSPECIAL(currentchar));
    gt_alphabet_echo_pretty_symbol(alpha,fpout,currentchar);
  }
}

void encseq2fastaoutput(FILE *fpout,
                        const char *desc,
                        const Encodedsequence *encseq,
                        Readmode readmode,
                        Seqpos start,
                        Seqpos wlen,
                        unsigned long width)
{
  gt_assert(width > 0);
  if (desc == NULL)
  {
    fprintf(fpout,">\n");
  } else
  {
    fprintf(fpout,">%s\n",desc);
  }
  encseq2symbolstring(fpout,
                      encseq,
                      readmode,
                      start,
                      wlen,
                      width);
}

int echodescriptionandsequence(const GtStrArray *filenametab,GtError *err)
{
  GtSeqIterator *seqit;
  char *desc = NULL;
  const GtUchar *sequence;
  unsigned long seqlen;
  bool haserr = false;
  int retval;

  seqit = gt_seqiterator_new(filenametab, err);
  if (!seqit)
    return -1;
  while (true)
  {
    retval = gt_seqiterator_next(seqit,
                              &sequence,
                              &seqlen,
                              &desc,
                              err);
    if (retval < 0)
    {
      haserr = true;
      break;
    }
    if (retval == 0)
    {
      break;
    }
    symbolstring2fasta(stdout,desc,NULL,sequence,seqlen,70UL);
    FREESPACE(desc);
  }
  gt_seqiterator_delete(seqit);
  return haserr ? -1 : 0;
}
