/*
  Copyright (C) 2007 David Ellinghaus <dellinghaus@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdio.h>

#include "libgtcore/env.h"

#include "libgtmatch/arraydef.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/encseq-def.h"
#include "libgtmatch/pos2seqnum.pr"

#include "repeattypes.h"
#include "ltrharvest-opt.h"

void showonstdout (
  char *s)
{
  printf ("# %s\n", s);
}

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
    // no repeats
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

int simpleexactselfmatchoutput (
  /*@unused@*/void *info,
  Seqpos len,
  Seqpos pos1,
  Seqpos pos2)
{
  Seqpos tmp;

  if (pos1 > pos2)
  {
    tmp = pos1;
    pos1 = pos2;
    pos2 = tmp;
  }

  printf(FormatSeqpos " " 
         FormatSeqpos " " 
	 FormatSeqpos "\n",
	 PRINTSeqposcast(len),
         PRINTSeqposcast(pos1),
         PRINTSeqposcast(pos2));

  return 0;
}

int simpleexactselfmatchstore (
  LTRharvestoptions *info,
  Seqpos len,
  Seqpos pos1,
  Seqpos pos2)
{
  Env *env = info->env;
  Seqpos tmp,
         totallength;
  const Encodedsequence *encseq = info->repeatinfo.suffixarrayptr->encseq;  
  unsigned long numofdbsequences = 
                        info->repeatinfo.suffixarrayptr->numofdbsequences;
  Seqpos *markpos;
  unsigned long i,
                contignumber = 0,
                seqnum1,
		seqnum2; 
  bool samecontig = false;

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
  // at least two db sequences
  else
  {
    markpos = calculatemarkpositions(encseq, numofdbsequences, env);
    if(markpos == NULL)
    {
      return -1;  
    }
    totallength = getencseqtotallength(encseq);
    for( i = 0; i < numofdbsequences - 1; i++)
    {
      seqnum1 = getrecordnum(markpos, numofdbsequences, 
	                totallength, pos1, env);
      if( seqnum1 == numofdbsequences)
      {
        return -1;
      }

      seqnum2 = getrecordnum(markpos, numofdbsequences, 
	                totallength, pos2, env);
      if( seqnum2 == numofdbsequences)
      {
        return -1;
      }

      if( seqnum1 == seqnum2 )
      {
/*
#ifdef DEBUG
	  printf("accepted:\n");
	  printf("pos1: " FormatSeqpos "\n", PRINTSeqposcast(pos1));
	  printf("pos2: " FormatSeqpos "\n", PRINTSeqposcast(pos2));
	  printf("i: " FormatSeqpos "\n", PRINTSeqposcast(i));
#endif
*/
	  samecontig = true;
          contignumber = seqnum1;
	  break;
      }
    }
    
    FREESPACE(markpos);
  }


  //test maximal length of candidate pair and distance constraints
  if( samecontig && (len <= (Seqpos) info->repeatinfo.lmax) && 
      ( (Seqpos) info->repeatinfo.dmin <= tmp) && 
        (tmp <= (Seqpos) info->repeatinfo.dmax) )
  {
    Repeat *nextfreerepeatptr;

    GETNEXTFREEINARRAY(nextfreerepeatptr, &info->repeatinfo.repeats, 
                       Repeat, 10);
#ifdef DEBUG
    printf("maximal repeat pos1: " FormatSeqpos "\n",
               PRINTSeqposcast(pos1));
    printf("maximal repeat pos2: " FormatSeqpos "\n",
               PRINTSeqposcast(pos2));
    printf("len: " FormatSeqpos "\n",
               PRINTSeqposcast(len));
    printf("seq number: %lu\n\n", contignumber);
#endif
    nextfreerepeatptr->pos1 = pos1;
    nextfreerepeatptr->offset = tmp;
    nextfreerepeatptr->len = len;
    nextfreerepeatptr->contignumber = contignumber;
  }

  return 0;
}

/*
int subsimpleexactselfmatchstore (
  SubRepeatInfo * info,
  Seqpos len,
  Seqpos pos1,
  Seqpos pos2)
{
  Seqpos tmp;

  if (pos1 > pos2)
  {
    tmp = pos1;
    pos1 = pos2;
    pos2 = tmp;
  }

  tmp = (pos2 - pos1);

  //test maximal length of candidate pair and if one instance per sequence
  if ((pos1 < (Uint) info->separatorpos)
      && ((Uint) info->separatorpos < pos2))
  {
    Repeat *nextfreerepeatptr;
#ifdef DEBUG
    printf("%lu %lu %lu\n",
	       PRINTSeqposcast(len),
	       PRINTSeqposcast(info->offset1 + pos1),
	       PRINTSeqposcast(info->offset2 - info->separatorpos - 1 + pos2));
#endif
    GETNEXTFREEINARRAY (nextfreerepeatptr, &info->repeats, Repeat, 10);
    nextfreerepeatptr->pos1 = pos1;
    nextfreerepeatptr->offset = tmp;
    nextfreerepeatptr->len = len;
  }

  return 0;
}
*/
