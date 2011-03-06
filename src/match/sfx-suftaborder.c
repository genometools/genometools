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
#include "sfx-suffixgetset.h"

static void showcomparisonfailure(const char *filename,
                                  int line,
                                  const char *where,
                                  const GtEncseq *encseq,
                                  GtReadmode readmode,
                                  const GtSuffixsortspace *suffixsortspace,
                                  unsigned long subbucketleft,
                                  unsigned long depth,
                                  unsigned long idx1,
                                  unsigned long idx2,
                                  int cmp,
                                  unsigned long maxlcp)
{
  unsigned long pos1, pos2;

  pos1 = gt_suffixsortspace_get(suffixsortspace,subbucketleft,idx1);
  pos2 = gt_suffixsortspace_get(suffixsortspace,subbucketleft,idx2);
  fprintf(stderr,"ERROR: file \"%s\", line %d: ",filename,line);
  fprintf(stderr,"%s(%lu vs %lu"
                 " %lu=\"",
                       where,
                       idx1,
                       idx2,
                       pos1);
  gt_encseq_showatstartposwithdepth(stderr,encseq,readmode,pos1,depth);
  fprintf(stderr,"\",\"");
  gt_encseq_showatstartposwithdepth(stderr,encseq,readmode,pos2,depth);
  fprintf(stderr,"\"=%lu)=%d with maxlcp %lu,depth=%lu\n",pos2,cmp,
          maxlcp,depth);
}

void gt_checkifprefixesareidentical(const char *filename,
                                    int line,
                                    const GtEncseq *encseq,
                                    GtReadmode readmode,
                                    const GtSuffixsortspace *suffixsortspace,
                                    unsigned long subbucketleft,
                                    unsigned long width,
                                    unsigned long depth)
{
  unsigned long idx, maxlcp, pos1, pos2;
  int cmp;
  GtEncseqReader *esr1, *esr2;

  esr1 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  esr2 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  gt_assert(depth > 0);
  for (idx = 0; idx < width-1; idx++)
  {
    pos1 = gt_suffixsortspace_get(suffixsortspace,subbucketleft,idx);
    pos2 = gt_suffixsortspace_get(suffixsortspace,subbucketleft,idx+1);
    cmp = gt_encseq_check_comparetwosuffixes(encseq,
                                             readmode,
                                             &maxlcp,
                                             false, /* specialsareequal */
                                             true,/* specialsareequalatdepth0 */
                                             depth,pos1,
                                             pos2,
                                             esr1,
                                             esr2);
    gt_assert(maxlcp <= depth);
    if (cmp != 0 || maxlcp < depth)
    {
      showcomparisonfailure(filename,
                            line,
                            "checkifprefixesareidentical",
                            encseq,
                            readmode,
                            suffixsortspace,
                            subbucketleft,
                            depth,
                            idx,idx+1,cmp,maxlcp);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
  gt_encseq_reader_delete(esr1);
  gt_encseq_reader_delete(esr2);
}

void gt_showentiresuftab(const GtEncseq *encseq,
                         GtReadmode readmode,
                         const GtSuffixsortspace *suffixsortspace,
                         unsigned long subbucketleft,
                         unsigned long depth)
{
  unsigned long idx, pos, totallength = gt_encseq_total_length(encseq);

  for (idx = 0; idx <= totallength; idx++)
  {
    pos = gt_suffixsortspace_get(suffixsortspace,subbucketleft,idx);
    printf("suftab[%lu]=%lu ",idx,pos);
    gt_encseq_showatstartposwithdepth(stdout,encseq,readmode,pos,depth);
    printf("\n");
  }
}

void gt_checksortedsuffixes(const char *filename,
                            int line,
                            const GtEncseq *encseq,
                            GtReadmode readmode,
                            const GtSuffixsortspace *suffixsortspace,
                            unsigned long subbucketleft,
                            unsigned long numberofsuffixes,
                            bool specialsareequal,
                            bool specialsareequalatdepth0,
                            unsigned long depth)
{
  unsigned long idx, pos1, pos2, maxlcp,
                totallength = gt_encseq_total_length(encseq);
  GtEncseqReader *esr1, *esr2;
  int cmp;

  gt_assert(!specialsareequal || specialsareequalatdepth0);
  esr1 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  esr2 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  gt_assert(numberofsuffixes > 0);
  pos1 = gt_suffixsortspace_get(suffixsortspace,subbucketleft,0);
  gt_assert(pos1 < totallength);
  for (idx = 1UL; idx < numberofsuffixes; idx++)
  {
    pos2 = gt_suffixsortspace_get(suffixsortspace,subbucketleft,idx);
    if (pos2 < totallength)
    {
      cmp = gt_encseq_check_comparetwosuffixes(encseq,
                                               readmode,
                                               &maxlcp,
                                               specialsareequal,
                                               specialsareequalatdepth0,
                                               depth,
                                               pos1,
                                               pos2,
                                               esr1,
                                               esr2);
      if (cmp > 0)
      {
        showcomparisonfailure(filename,
                              line,
                              "checksortedsuffixes",
                              encseq,
                              readmode,
                              suffixsortspace,
                              subbucketleft,
                              depth,
                              idx-1,
                              idx,
                              cmp,
                              maxlcp);
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
      gt_assert(depth == 0 || maxlcp <= depth);
    }
    pos1 = pos2;
  }
  gt_encseq_reader_delete(esr1);
  gt_encseq_reader_delete(esr2);
}
