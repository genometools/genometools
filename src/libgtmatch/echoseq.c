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

#include <assert.h>
#include "libgtcore/chardef.h"
#include "libgtcore/error.h"
#include "libgtcore/seqiterator.h"
#include "spacedef.h"
#include "alphadef.h"
#include "encseq-def.h"

unsigned long *calcdescendpositions(const char *destab,
                                    unsigned long destablength,
                                    unsigned long numofsequences)
{
  unsigned long *descendtab, i, idx = 0;

  ALLOCASSIGNSPACE(descendtab,NULL,unsigned long,numofsequences);
  assert(destab != NULL);
  for (i=0; i<destablength; i++)
  {
    if (destab[i] == '\n')
    {
      assert(idx < numofsequences);
      descendtab[idx++] = i;
    }
  }
  assert(idx == numofsequences);
  return descendtab;
}

const char *retriesequencedescription(unsigned long *desclen,
                                      const char *destab,
                                      const unsigned long *descendtab,
                                      unsigned long seqnum)
{
  if (seqnum == 0)
  {
    *desclen = descendtab[0];
    return destab;
  }
  assert(descendtab[seqnum-1] < descendtab[seqnum]);
  *desclen = descendtab[seqnum] - descendtab[seqnum-1] - 1;
  return destab + descendtab[seqnum-1] + 1;
}

void checkalldescriptions(const char *destab,unsigned long destablength,
                          unsigned long numofsequences)
{
  unsigned long *descendtab, desclen, seqnum, totaldesclength, offset = 0;
  const char *desptr;
  char *copydestab;

  descendtab = calcdescendpositions(destab,
                                    destablength,
                                    numofsequences);
  totaldesclength = numofsequences; /* for each new line */
  for (seqnum = 0; seqnum < numofsequences; seqnum++)
  {
    desptr = retriesequencedescription(&desclen,
                                       destab,
                                       descendtab,
                                       seqnum);
    totaldesclength += desclen;
  }
  ALLOCASSIGNSPACE(copydestab,NULL,char,totaldesclength);
  for (seqnum = 0; seqnum < numofsequences; seqnum++)
  {
    desptr = retriesequencedescription(&desclen,
                                       destab,
                                       descendtab,
                                       seqnum);
    strncpy(copydestab + offset,desptr,(size_t) desclen);
    copydestab[offset+desclen] = '\n';
    offset += (desclen+1);
  }
  if (strncmp(copydestab,destab,(size_t) totaldesclength) != 0)
  {
    fprintf(stderr,"different descriptions\n");
    exit(EXIT_FAILURE); /* Programm error */
  }
  FREESPACE(copydestab);
  FREESPACE(descendtab);
}

void symbolstring2fasta(FILE *fpout,
                        const char *desc,
                        const Alphabet *alpha,
                        const Uchar *w,
                        unsigned long wlen,
                        unsigned long width)
{
  unsigned long i, j;
  Uchar currentchar;

  assert(width > 0);
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
    if (currentchar == (Uchar) SEPARATOR)
    {
      fprintf(fpout,"\n>\n");
      j = 0;
    } else
    {
      echoprettysymbol(fpout,alpha,currentchar);
    }
    if (i == wlen - 1)
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

void encseq2symbolstring(FILE *fpout,
                         const Alphabet *alpha,
                         const Encodedsequence *encseq,
                         Readmode readmode,
                         Seqpos start,
                         unsigned long wlen,
                         unsigned long width)
{
  unsigned long j;
  Seqpos idx, lastpos;
  Uchar currentchar;
  Encodedsequencescanstate *esr;

  esr = newEncodedsequencescanstate();
  initEncodedsequencescanstate(esr,encseq,readmode,start);
  assert(width > 0);
  lastpos = start + wlen - 1;
  for (idx = start, j = 0; /* Nothing */ ; idx++)
  {
    currentchar = sequentialgetencodedchar(encseq,esr,idx,readmode);
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
  freeEncodedsequencescanstate(&esr);
}

void encseq2fastaoutput(FILE *fpout,
                        const char *desc,
                        const Alphabet *alpha,
                        const Encodedsequence *encseq,
                        Readmode readmode,
                        Seqpos start,
                        unsigned long wlen,
                        unsigned long width)
{
  assert(width > 0);
  if (desc == NULL)
  {
    fprintf(fpout,">\n");
  } else
  {
    fprintf(fpout,">%s\n",desc);
  }
  encseq2symbolstring(fpout,
                      alpha,
                      encseq,
                      readmode,
                      start,
                      wlen,
                      width);
}

int echodescriptionandsequence(const StrArray *filenametab,Error *err)
{
  SeqIterator *seqit;
  char *desc = NULL;
  const Uchar *sequence;
  unsigned long seqlen;
  bool haserr = false;
  int retval;

  seqit = seqiterator_new(filenametab,NULL,true);
  while (true)
  {
    retval = seqiterator_next(seqit,
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
  seqiterator_delete(seqit);
  return haserr ? -1 : 0;
}
