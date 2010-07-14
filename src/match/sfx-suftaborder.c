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

#include <stdio.h>
#include <stdlib.h>
#include "core/unused_api.h"
#include "core/chardef.h"
#include "core/minmax.h"
#include "core/encseq.h"
#include "sfx-suftaborder.h"
#include "suffixptr.h"

static void showcomparisonfailure(const char *filename,
                                  int line,
                                  const char *where,
                                  const GtEncseq *encseq,
                                  GtReadmode readmode,
                                  const Suffixptr *suftab,
                                  unsigned long depth,
                                  unsigned long idx1,
                                  unsigned long idx2,
                                  int cmp,
                                  unsigned long maxlcp)
{
  fprintf(stderr,"ERROR: file \"%s\", line %d: ",filename,line);
  fprintf(stderr,"%s(%lu vs %lu"
                 " %lu=\"",
                       where,
                       idx1,
                       idx2,
                       SUFFIXPTRGET(suftab,idx1));
  gt_encseq_showatstartposwithdepth(stderr,encseq,readmode,
                                    SUFFIXPTRGET(suftab,idx1),depth);
  fprintf(stderr,"\",\"");
  gt_encseq_showatstartposwithdepth(stderr,encseq,readmode,
                                    SUFFIXPTRGET(suftab,idx2),depth);
  fprintf(stderr,"\"=%lu)=%d with maxlcp %lu\n",SUFFIXPTRGET(suftab,idx2),
                                                cmp,
                                                maxlcp);
}

void gt_checkifprefixesareidentical(const char *filename,
                                    int line,
                                    const GtEncseq *encseq,
                                    GtReadmode readmode,
                                    const Suffixptr *suftab,
                                    unsigned int prefixlength,
                                    unsigned long depth,
                                    unsigned long left,
                                    unsigned long right)
{
  unsigned long idx, maxlcp;
  int cmp;
  GtEncseqReader *esr1, *esr2;
  bool haserr = false;

  esr1 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  esr2 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  for (idx = left; idx < right; idx++)
  {
    cmp = gt_encseq_comparetwosuffixes(encseq,
                                       readmode,
                                       &maxlcp,
                                       false,
                                       true,
                                       depth,
                                       SUFFIXPTRGET(suftab,idx),
                                       SUFFIXPTRGET(suftab,idx+1),
                                       esr1,
                                       esr2);
    if (cmp != 0 || maxlcp != (unsigned long) prefixlength)
    {
      showcomparisonfailure(filename,
                            line,
                            "checkifprefixesareidentical",
                            encseq,
                            readmode,
                            suftab,
                            depth,
                            idx,idx+1,cmp,maxlcp);
      haserr = true;
      break;
    }
  }
  gt_encseq_reader_delete(esr1);
  gt_encseq_reader_delete(esr2);
  if (haserr)
  {
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

void gt_showentiresuftab(const GtEncseq *encseq,
                         GtReadmode readmode,
                         const Suffixptr *suftab,
                         unsigned long depth)
{
  unsigned long idx, totallength = gt_encseq_total_length(encseq);

  for (idx = 0; idx <= totallength; idx++)
  {
    printf("suftab[%lu]=%lu ",
            idx,SUFFIXPTRGET(suftab,idx));
    gt_encseq_showatstartposwithdepth(stdout,encseq,readmode,
                                      SUFFIXPTRGET(suftab,idx),depth);
    printf("\n");
  }
}

void gt_checksortedsuffixes(const char *filename,
                            int line,
                            const GtEncseq *encseq,
                            GtReadmode readmode,
                            const Suffixptr *subbucket,
                            GT_UNUSED unsigned long subbucketleft,
                            unsigned long numberofsuffixes,
                            bool specialsareequal,
                            bool specialsareequalatdepth0,
                            unsigned long depth)
{
  unsigned long idx, maxlcp, totallength = gt_encseq_total_length(encseq);
  GtEncseqReader *esr1, *esr2;
  int cmp;

  gt_assert(!specialsareequal || specialsareequalatdepth0);
  esr1 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  esr2 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  gt_assert(numberofsuffixes > 0);
  gt_assert(SUFFIXPTRGET(subbucket,0) < totallength);
  for (idx = 1UL; idx < numberofsuffixes; idx++)
  {
    if (idx < numberofsuffixes - 1)
    {
      gt_assert(SUFFIXPTRGET(subbucket,idx) < totallength);
      cmp = gt_encseq_comparetwosuffixes(encseq,
                                         readmode,
                                         &maxlcp,
                                         specialsareequal,
                                         specialsareequalatdepth0,
                                         depth,
                                         SUFFIXPTRGET(subbucket,idx-1),
                                         SUFFIXPTRGET(subbucket,idx),
                                         esr1,
                                         esr2);
      if (cmp > 0)
      {
        showcomparisonfailure(filename,
                              line,
                              "checksortedsuffixes",
                              encseq,
                              readmode,
                              subbucket,
                              depth,
                              idx-1,
                              idx,
                              cmp,
                              maxlcp);
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
    } else
    {
      if (numberofsuffixes == totallength+1)
      {
        gt_assert(SUFFIXPTRGET(subbucket,idx) == totallength);
      }
    }
  }
  gt_encseq_reader_delete(esr1);
  gt_encseq_reader_delete(esr2);
}
