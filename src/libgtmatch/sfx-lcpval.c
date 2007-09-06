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

#include "seqpos-def.h"
#include "encseq-def.h"
#include "spacedef.h"
#include "sfx-lcpval.h"

#include "sfx-cmpsuf.pr"

 struct Lcpvalueiterator
{
  Seqpos relpos, lastsuftabentryofpreviouspart;
  Readmode readmode;
  const Encodedsequence *encseq;
};

Lcpvalueiterator *newLcpvalueiterator(const Encodedsequence *encseq,
                                      Readmode readmode,
                                      Env *env)
{
  Lcpvalueiterator *lvi;
  
  ALLOCASSIGNSPACE(lvi,NULL,Lcpvalueiterator,1);
  lvi->encseq = encseq;
  lvi->relpos = 0;
  lvi->readmode = readmode;
  lvi->lastsuftabentryofpreviouspart = 0;
  return lvi;
}

Seqpos nextLcpvalueiterator(Lcpvalueiterator *lvi,
                            bool firstpage,
                            const Seqpos *suftabptr,
                            Seqpos numberofsuffixes)
{
  Seqpos lcpvalue;

  assert(lvi->relpos < numberofsuffixes);
  if(firstpage && lvi->relpos == 0)
  {
    lcpvalue = 0;
  } else
  {
    int cmp;

    cmp = comparetwosuffixes(lvi->encseq,
                             lvi->readmode,
                             &lcpvalue,
                             false,
                             false,
                             0,
                             lvi->relpos > 0 
                               ? suftabptr[lvi->relpos-1]
                               : lvi->lastsuftabentryofpreviouspart,
                             suftabptr[lvi->relpos]);
    if (cmp > 0)
    {
      fprintf(stderr,"pos = " FormatSeqpos
              ": cmp " FormatSeqpos
              " " FormatSeqpos " = %d",
              PRINTSeqposcast(lvi->relpos),
              lvi->relpos > 0 
                ? PRINTSeqposcast(suftabptr[lvi->relpos-1])
                : PRINTSeqposcast(lvi->lastsuftabentryofpreviouspart),
              PRINTSeqposcast(suftabptr[lvi->relpos]),
              cmp);
      exit(EXIT_FAILURE);
    }
  }
  if(lvi->relpos + 1 == numberofsuffixes)
  {
    lvi->relpos = 0;
    lvi->lastsuftabentryofpreviouspart = suftabptr[numberofsuffixes-1];
  } else
  {
    lvi->relpos++;
  }
  return lcpvalue;
}

void freeLcpvalueiterator(Lcpvalueiterator **lvi,Env *env)
{
  FREESPACE(*lvi);
}
