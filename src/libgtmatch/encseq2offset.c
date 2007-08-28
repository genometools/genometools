/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "spacedef.h"
#include "chardef.h"
#include "sarr-def.h"
#include "safecast-gen.h"

 DECLARESAFECASTFUNCTION(uint64_t,uint64_t,Seqpos,Seqpos)

Seqpos *encseqtable2seqoffsets(Seqpos *totallength,
                               Specialcharinfo *specialcharinfo,
                               const Suffixarray *suffixarraytable,
                               uint32_t numofindexes,
                               Env *env)
{
  uint32_t idx;
  Uchar lastofprevious, firstofcurrent;
  Seqpos tmplength, *sequenceoffsettable;
  uint64_t tmpspecialcharacters,
           tmpspecialranges,
           tmplarge;

  env_error_check(env);
  assert(numofindexes > 0);
  ALLOCASSIGNSPACE(sequenceoffsettable,NULL,Seqpos,numofindexes);
  tmpspecialcharacters = (uint64_t) (numofindexes-1);
  tmpspecialranges = 0;
  for (idx=0; idx<numofindexes; idx++)
  {
    if (idx == 0)
    {
      tmplength = 0;
      sequenceoffsettable[idx] = 0;
    } else
    {
      tmplength = getencseqtotallength(suffixarraytable[idx - 1].encseq);
      sequenceoffsettable[idx]
        = sequenceoffsettable[idx-1] + tmplength + (Seqpos) 1;
    }
    tmpspecialcharacters
      += (uint64_t) suffixarraytable[idx].specialcharinfo.specialcharacters;
    tmpspecialranges
      += (uint64_t) suffixarraytable[idx].specialcharinfo.specialranges;
    if (idx > 0)
    {
      lastofprevious = getencodedchar(suffixarraytable[idx - 1].encseq,
                                      tmplength-1,
                                      suffixarraytable[idx - 1].readmode);
      firstofcurrent = getencodedchar(suffixarraytable[idx].encseq,
                                      0,
                                      suffixarraytable[idx].readmode);
      if (ISSPECIAL(lastofprevious))
      {
         if (ISSPECIAL(firstofcurrent))
         {
           tmpspecialranges--;
         }
      } else
      {
        if (ISNOTSPECIAL(firstofcurrent))
        {
          tmpspecialranges++;
        }
      }
    }
    tmplarge = (uint64_t) sequenceoffsettable[idx] +
               (uint64_t) getencseqtotallength(suffixarraytable[idx].encseq);
    (void) CALLCASTFUNC(uint64_t,Seqpos,tmplarge);
    (void) CALLCASTFUNC(uint64_t,Seqpos,tmpspecialcharacters);
    (void) CALLCASTFUNC(uint64_t,Seqpos,tmpspecialranges);
    printf("# seqlen[%u] = " FormatSeqpos "\n",
           (unsigned int) idx,
           PRINTSeqposcast(getencseqtotallength(suffixarraytable[idx].encseq)));
  }
  tmplength = getencseqtotallength(suffixarraytable[numofindexes -1].encseq);
  *totallength = sequenceoffsettable[numofindexes-1] + tmplength;
  specialcharinfo->specialcharacters = (Seqpos) tmpspecialcharacters;
  specialcharinfo->specialranges = (Seqpos) tmpspecialranges;
  specialcharinfo->lengthofspecialprefix
    = suffixarraytable[0].specialcharinfo.lengthofspecialprefix;
  specialcharinfo->lengthofspecialsuffix
    = suffixarraytable[idx-1].specialcharinfo.lengthofspecialsuffix;
  return sequenceoffsettable;
}
