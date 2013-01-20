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

#include "core/encseq.h"
#include "core/unused_api.h"
#include "esa-lcpval.h"

 struct Lcpvalueiterator
{
  unsigned long relpos,
         lastsuftabentry;
  GtReadmode readmode;
  const GtEncseq *encseq;
  GtEncseqReader *esr1, *esr2;
};

Lcpvalueiterator *gt_newLcpvalueiterator(const GtEncseq *encseq,
                                         GtReadmode readmode)
{
  Lcpvalueiterator *lvi;

  lvi = gt_malloc(sizeof *lvi);
  lvi->esr1 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  lvi->esr2 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  lvi->encseq = encseq;
  lvi->relpos = 0;
  lvi->readmode = readmode;
  lvi->lastsuftabentry = 0;
  return lvi;
}

unsigned long gt_nextLcpvalueiterator(Lcpvalueiterator *lvi,
                                      bool firstpage,
                                      const ESASuffixptr *suftabptr,
                                      unsigned long numberofsuffixes)
{
  unsigned long lcpvalue;

  gt_assert(lvi->relpos < numberofsuffixes);
  if (firstpage && lvi->relpos == 0)
  {
    lcpvalue = 0;
  } else
  {
    GT_UNUSED int cmp;

    cmp = gt_encseq_check_comparetwosuffixes(lvi->encseq,
                                             lvi->readmode,
                                             &lcpvalue,
                                             false,
                                             false,
                                             0,
                                             lvi->lastsuftabentry,
                                             ESASUFFIXPTRGET(suftabptr,
                                                             lvi->relpos),
                                             lvi->esr1,
                                             lvi->esr2);
#ifndef NDEBUG
    if (cmp > 0)
    {
      fprintf(stderr,"pos=%lu"
              ": cmp %lu"
              " %lu = %d, lcpval=%lu\n",
              lvi->relpos,
              lvi->lastsuftabentry,
              ESASUFFIXPTRGET(suftabptr,lvi->relpos),
              cmp,
              lcpvalue);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
#endif
  }
  lvi->lastsuftabentry = ESASUFFIXPTRGET(suftabptr,lvi->relpos);
  if (lvi->relpos + 1 == numberofsuffixes)
  {
    lvi->relpos = 0;
  } else
  {
    lvi->relpos++;
  }
  return lcpvalue;
}

void gt_freeLcpvalueiterator(Lcpvalueiterator *lvi)
{
  if (lvi != NULL)
  {
    gt_encseq_reader_delete(lvi->esr1);
    gt_encseq_reader_delete(lvi->esr2);
    gt_free(lvi);
  }
}
