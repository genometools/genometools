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
#include "spacedef.h"
#include "sarr-def.h"
#include "core/safecast-gen.h"

unsigned long *encseqtable2sequenceoffsets(unsigned long *totallength,
                                    GtSpecialcharinfo *specialcharinfo,
                                    const Suffixarray *suffixarraytable,
                                    unsigned int numofindexes)
{
  unsigned int idx;
  GtUchar lastofprevious, firstofcurrent;
  unsigned long tmplength, *sequenceoffsettable;
  uint64_t tmpspecialcharacters,
           tmpspecialranges,
           tmprealspecialranges,
           tmplarge;

  gt_assert(numofindexes > 0);
  ALLOCASSIGNSPACE(sequenceoffsettable,NULL,unsigned long,numofindexes);
  tmpspecialcharacters = (uint64_t) (numofindexes-1);
  tmpspecialranges = 0;
  tmprealspecialranges = 0;
  for (idx=0; idx<numofindexes; idx++)
  {
    if (idx == 0)
    {
      tmplength = 0;
      sequenceoffsettable[idx] = 0;
    } else
    {
      tmplength =
              gt_encodedsequence_total_length(suffixarraytable[idx - 1].encseq);
      sequenceoffsettable[idx]
        = sequenceoffsettable[idx-1] + tmplength + (unsigned long) 1;
    }
    tmpspecialcharacters
      += (uint64_t) gt_encodedsequence_specialcharacters(
                                                  suffixarraytable[idx].encseq);
    tmpspecialranges
      += (uint64_t) gt_encodedsequence_specialranges(
                                                  suffixarraytable[idx].encseq);
    tmprealspecialranges
      += (uint64_t) gt_encodedsequence_realspecialranges(
                                                  suffixarraytable[idx].encseq);
    if (idx > 0)
    {
      /* Random access */
      lastofprevious
        = gt_encodedsequence_getencodedchar(suffixarraytable[idx - 1].encseq,
                                            tmplength-1,
                                            suffixarraytable[idx - 1].readmode);
      /* Random access */
      firstofcurrent
        = gt_encodedsequence_getencodedchar(suffixarraytable[idx].encseq,
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
    tmplarge = (uint64_t) sequenceoffsettable[idx] +
       (uint64_t) gt_encodedsequence_total_length(suffixarraytable[idx].encseq);
    (void) CALLCASTFUNC(uint64_t,unsigned_long,tmplarge);
    (void) CALLCASTFUNC(uint64_t,unsigned_long,tmpspecialcharacters);
    (void) CALLCASTFUNC(uint64_t,unsigned_long,tmpspecialranges);
    (void) CALLCASTFUNC(uint64_t,unsigned_long,tmprealspecialranges);
    printf("# seqlen[%u] = %lu\n",
           idx,
           gt_encodedsequence_total_length(suffixarraytable[idx].encseq));
  }
  tmplength
    = gt_encodedsequence_total_length(suffixarraytable[numofindexes -1].encseq);
  *totallength = sequenceoffsettable[numofindexes-1] + tmplength;
  specialcharinfo->specialcharacters = (unsigned long) tmpspecialcharacters;
  specialcharinfo->specialranges = (unsigned long) tmpspecialranges;
  specialcharinfo->realspecialranges = (unsigned long) tmprealspecialranges;
  specialcharinfo->lengthofspecialprefix
    = gt_encodedsequence_lengthofspecialprefix(suffixarraytable[0].encseq);
  specialcharinfo->lengthofspecialsuffix
    = gt_encodedsequence_lengthofspecialsuffix(suffixarraytable[idx-1].encseq);
  return sequenceoffsettable;
}
