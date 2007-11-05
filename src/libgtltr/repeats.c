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

#include <stdio.h>
#include "libgtcore/env.h"
#include "libgtcore/arraydef.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/encseq-def.h"
#include "libgtmatch/intcode-def.h"
#include "repeattypes.h"
#include "ltrharvest-opt.h"

#include "libgtmatch/pos2seqnum.pr"

void showrepeats (
  RepeatInfo * repeatinfo,
  unsigned long seedminlength)
{
  ArrayRepeat *repeats = &repeatinfo->repeats;
  Repeat *reptab = repeats->spaceRepeat;
  unsigned long i;

  printf("# there are %lu seed-pairs as candidates for LTRs\n"
         "# of user defined minimal seedlength %lu:",
          repeats->nextfreeRepeat, seedminlength);
  if (repeats->nextfreeRepeat == (unsigned long)0)
  {
    /* no repeats */
    printf ("\nnone");
  }
  else
  {
    for (i = 0; i < repeats->nextfreeRepeat; i++)
    {
      printf ("\n#   len: " FormatSeqpos " ", PRINTSeqposcast(reptab[i].len));
      printf ("pos1: " FormatSeqpos " ", PRINTSeqposcast(reptab[i].pos1));
      printf ("pos2: " FormatSeqpos " ",
              PRINTSeqposcast((reptab[i].pos1 + reptab[i].offset)));
      printf ("offset: " FormatSeqpos " ", PRINTSeqposcast(reptab[i].offset));
    }
  }
  printf ("\n");
}

int simpleexactselfmatchstore (
  LTRharvestoptions *info,
  Seqpos len,
  Seqpos pos1,
  Seqpos pos2,
  Env *env)
{
  Seqpos tmp,
         totallength;
  const Encodedsequence *encseq =
          encseqSequentialsuffixarrayreader(info->repeatinfo.ssarptr);
  unsigned long numofdbsequences =
       numofdbsequencesSequentialsuffixarrayreader(info->repeatinfo.ssarptr);
  Seqpos *markpos = info->markpos;
  unsigned long i,
                contignumber = 0,
                seqnum1,
		seqnum2;
  bool samecontig = false;

  env_error_check(env);

  if (pos1 > pos2)
  {
    tmp = pos1;
    pos1 = pos2;
    pos2 = tmp;
  }

  tmp = (pos2 - pos1);
  if ( numofdbsequences < (unsigned long) 2 )
  {
    samecontig = true;
    contignumber = (unsigned long)0;
  }
  /* at least two db sequences */
  else
  {
    if (markpos == NULL)
    {
      return -1;
    }
    totallength = getencseqtotallength(encseq);
    for ( i = 0; i < numofdbsequences - 1; i++)
    {
      seqnum1 = getrecordnumSeqpos(markpos, numofdbsequences,
	                totallength, pos1, env);
      if ( seqnum1 == numofdbsequences)
      {
        return -1;
      }

      seqnum2 = getrecordnumSeqpos(markpos, numofdbsequences,
	                totallength, pos2, env);
      if ( seqnum2 == numofdbsequences)
      {
        return -1;
      }

      if ( seqnum1 == seqnum2 )
      {
	  env_log_log(env, "accepted:\n");
	  env_log_log(env, "pos1: " FormatSeqpos "\n", PRINTSeqposcast(pos1));
	  env_log_log(env, "pos2: " FormatSeqpos "\n", PRINTSeqposcast(pos2));
	  env_log_log(env, "i: %lu\n", i);
	  samecontig = true;
          contignumber = seqnum1;
	  break;
      }
    }
  }

  /*test maximal length of candidate pair and distance constraints*/
  if ( samecontig && (len <= (Seqpos) info->repeatinfo.lmax) &&
      ( (Seqpos) info->repeatinfo.dmin <= tmp) &&
        (tmp <= (Seqpos) info->repeatinfo.dmax) )
  {
    Repeat *nextfreerepeatptr;

    GETNEXTFREEINARRAY(nextfreerepeatptr, &info->repeatinfo.repeats,
                       Repeat, 10);
    env_log_log(env, "maximal repeat pos1: " FormatSeqpos "\n",
               PRINTSeqposcast(pos1));
    env_log_log(env, "maximal repeat pos2: " FormatSeqpos "\n",
               PRINTSeqposcast(pos2));
    env_log_log(env, "len: " FormatSeqpos "\n",
               PRINTSeqposcast(len));
    env_log_log(env, "seq number: %lu\n\n", contignumber);
    nextfreerepeatptr->pos1 = pos1;
    nextfreerepeatptr->offset = tmp;
    nextfreerepeatptr->len = len;
    nextfreerepeatptr->contignumber = contignumber;
  }

  return 0;
}

int subsimpleexactselfmatchstore (
  void *info,
  unsigned long len,
  Seqpos dbstart,
   /*@unused@*/ uint64_t queryoffset,
  unsigned long querystart,
  /*@unused@*/ Env *env)
{
  Repeat *nextfreerepeatptr;
  SubRepeatInfo *sri = (SubRepeatInfo *) info;

  GETNEXTFREEINARRAY (nextfreerepeatptr, &sri->repeats, Repeat, 10);
  nextfreerepeatptr->pos1 = sri->offset1 + dbstart;
  nextfreerepeatptr->offset = sri->offset2 + (Seqpos)querystart -
    (sri->offset1 + dbstart);
  nextfreerepeatptr->len = (Seqpos)len;

  return 0;
}
