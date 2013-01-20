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

#include "core/chardef.h"
#include "core/ma_api.h"
#include "sarr-def.h"
#include "core/safecast-gen.h"

unsigned long *gt_encseqtable2sequenceoffsets(
                                    unsigned long *totallength,
                                    GtSpecialcharinfo *specialcharinfo,
                                    const Suffixarray *suffixarraytable,
                                    unsigned int numofindexes)
{
  unsigned int idx;
  GtUchar lastofprevious, firstofcurrent;
  unsigned long tmplength, numofsequences = 0, *sequenceoffsettable;
  uint64_t tmpspecialcharacters,
           tmpwildcards,
           tmpspecialranges,
           tmpwildcardranges,
           tmprealspecialranges,
           tmprealwildcardranges,
           tmpoffset;

  gt_assert(numofindexes > 0);
  sequenceoffsettable = gt_malloc(sizeof (*sequenceoffsettable) * numofindexes);
  tmpspecialcharacters = (uint64_t) (numofindexes-1);
  tmpwildcards = 0;
  tmpspecialranges = 0;
  tmpwildcardranges = 0;
  tmprealspecialranges = 0;
  tmprealwildcardranges = 0;
  for (idx=0; idx<numofindexes; idx++)
  {
    if (idx == 0)
    {
      tmplength = 0;
      sequenceoffsettable[idx] = 0;
    } else
    {
      tmplength = gt_encseq_total_length(suffixarraytable[idx - 1].encseq);
      sequenceoffsettable[idx] = sequenceoffsettable[idx-1] + tmplength + 1UL;
    }
    numofsequences += gt_encseq_num_of_sequences(suffixarraytable[idx].encseq);
    tmpspecialcharacters
      += (uint64_t) gt_encseq_specialcharacters(suffixarraytable[idx].encseq);
    tmpwildcards
      += (uint64_t) gt_encseq_wildcards(suffixarraytable[idx].encseq);
    tmpspecialranges
      += (uint64_t) gt_encseq_specialranges(suffixarraytable[idx].encseq);
    tmpwildcardranges
      += (uint64_t) gt_encseq_wildcardranges(suffixarraytable[idx].encseq);
    tmprealspecialranges
      += (uint64_t) gt_encseq_realspecialranges(suffixarraytable[idx].encseq);
    tmprealwildcardranges
      += (uint64_t) gt_encseq_realwildcardranges(suffixarraytable[idx].encseq);
    if (idx > 0)
    {
      /* Random access */
      lastofprevious
        = gt_encseq_get_encoded_char(suffixarraytable[idx-1].encseq,
                                     tmplength-1,
                                     suffixarraytable[idx - 1].readmode);
      /* Random access */
      firstofcurrent
        = gt_encseq_get_encoded_char(suffixarraytable[idx].encseq,
                                     0,
                                     suffixarraytable[idx].readmode);
      if (ISSPECIAL(lastofprevious))
      {
        if (ISSPECIAL(firstofcurrent))
        {
          tmpspecialranges--;
          tmprealspecialranges--;
        }
      } else
      {
        if (ISNOTSPECIAL(firstofcurrent))
        {
          tmpspecialranges++;
          tmprealspecialranges++;
        }
      }
    }
    tmpoffset = (uint64_t) sequenceoffsettable[idx] +
       (uint64_t) gt_encseq_total_length(suffixarraytable[idx].encseq);
    (void) CALLCASTFUNC(uint64_t,unsigned_long,tmpoffset);
    (void) CALLCASTFUNC(uint64_t,unsigned_long,tmpspecialcharacters);
    (void) CALLCASTFUNC(uint64_t,unsigned_long,tmpwildcards);
    (void) CALLCASTFUNC(uint64_t,unsigned_long,tmpspecialranges);
    (void) CALLCASTFUNC(uint64_t,unsigned_long,tmpwildcardranges);
    (void) CALLCASTFUNC(uint64_t,unsigned_long,tmprealspecialranges);
    (void) CALLCASTFUNC(uint64_t,unsigned_long,tmprealwildcardranges);
    printf("# seqlen[%u] = %lu\n",
           idx,
           gt_encseq_total_length(suffixarraytable[idx].encseq));
  }
  tmplength = gt_encseq_total_length(suffixarraytable[numofindexes -1].encseq);
  *totallength = sequenceoffsettable[numofindexes-1] + tmplength;
  specialcharinfo->specialcharacters = (unsigned long) tmpspecialcharacters;
  specialcharinfo->wildcards = (unsigned long) tmpwildcards;
  specialcharinfo->specialranges = (unsigned long) tmpspecialranges;
  specialcharinfo->wildcardranges = (unsigned long) tmpwildcardranges;
  specialcharinfo->realspecialranges = (unsigned long) tmprealspecialranges;
  specialcharinfo->realwildcardranges = (unsigned long) tmprealwildcardranges;
  specialcharinfo->lengthofspecialprefix
    = gt_encseq_lengthofspecialprefix(suffixarraytable[0].encseq);
  specialcharinfo->lengthofwildcardprefix
    = gt_encseq_lengthofwildcardprefix(suffixarraytable[0].encseq);
  specialcharinfo->lengthofspecialsuffix
    = gt_encseq_lengthofspecialsuffix(suffixarraytable[idx-1].encseq);
  specialcharinfo->lengthofwildcardsuffix
    = gt_encseq_lengthofwildcardsuffix(suffixarraytable[idx-1].encseq);
  gt_assert(numofsequences > 0);
#ifndef NDEBUG
  gt_GtSpecialcharinfo_check(specialcharinfo,numofsequences - 1);
#endif
  return sequenceoffsettable;
}
