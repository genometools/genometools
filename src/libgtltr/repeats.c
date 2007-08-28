/* Prof.Kurtz files */
#include "libgtmatch/arraydef.h"
//#include "spacedef.h"

#include "libgtmatch/sarr-def.h"
#include "libgtmatch/encseq-def.h"

/* my files */
#include "repeattypes.h"

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
/*
int simpleexactselfmatchstore (
  RepeatInfo *info,
  Seqpos len,
  Seqpos pos1,
  Seqpos pos2)
{
  Seqpos tmp;
  unsigned int long;
  bool samecontig = false;
  const Encodedsequence *encseq = lo->repeatinfo.suffixarrayptr->encseq;  
  unsigned long numofdbsequences = lo->repeatinfo.suffixarrayptr->numofdbsequences;
  Seqpos *recordseps;
  //Multiseq *multiseqptr = &(info->virtualtreeptr->multiseq);
  //unsigned int numofsequences = multiseqptr->numofsequences;
  unsigned long contignumber = 0; 

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
    contignumber = 0;
  }
  // at least two db sequences
  else
  {
    recordseps = calculatemarkpositions(encseq, numofdbsequences, env);
    if(recordseps == NULL)
    {
      return -1;  
    }
    getrecordnum(recordseps, numofdbsequences-1, 
                 encseq->totallength, pos1, env);
    for( i = 0; i < numofdbsequences - 1; i++)
    {
      if( getrecordnum(recordseps, numofdbsequences-1, 
	                encseq->totallength, pos1, env) !=
	  getrecordnum(recordseps, numofdbsequences-1, 
	                encseq->totallength, pos2, env))
*/      /*if( (i == 0 || multiseqptr->markpos.spaceUint[i-1] < pos1)
	  && pos1 < multiseqptr->markpos.spaceUint[i])
      {
	if( (i == 0 || pos2 > multiseqptr->markpos.spaceUint[i-1])
	    && pos2 < multiseqptr->markpos.spaceUint[i])
	{*/
/*      {
#ifdef DEBUG
	  printf("accepted:\n");
	  printf("pos1: %lu\n", PRINTSeqposcast(pos1));
	  printf("pos2: %lu\n", PRINTSeqposcast(pos2));
	  printf("i: %lu\n", PRINTSeqposcast(i));
	  if(i > 0)
	  {
	    printf("separator pos at i-1: %ld\n",
	            PRINTSeqposcast(multiseqptr->markpos.spaceUint[i-1]));
	  }
	  printf("separator pos at i: %ld\n",
	            PRINTSeqposcast(multiseqptr->markpos.spaceUint[i]));
#endif
	  samecontig = true;
          contignumber = i;
	  break;
      }
    }
    if(pos1 > multiseqptr->markpos.spaceUint[numofsequences - 2])
    {
      if( pos2 > multiseqptr->markpos.spaceUint[numofsequences - 2] ) 
      {
	samecontig = True;
	contignumber = numofsequences - 1;
      }
    }
  }

  //test maximal length of candidate pair and distance constraints
  if (samecontig && (len <= info->lmax) && (info->dmin <= tmp)
      && (tmp <= info->dmax))
  {
    Repeat *nextfreerepeatptr;

    GETNEXTFREEINARRAY (nextfreerepeatptr, &info->repeats, Repeat, 10);
#ifdef DEBUG
    fprintf("maximal repeat pos1: %lu\n",
               PRINTSeqposcast(pos1));
    fprintf("maximal repeat pos2: %lu\n",
               PRINTSeqposcast(pos1 + tmp));
    fprintf("seq number: %lu\n\n",
               (Showuint) contignumber);
#endif
    nextfreerepeatptr->pos1 = pos1;
    nextfreerepeatptr->offset = tmp;
    nextfreerepeatptr->len = len;
    nextfreerepeatptr->contignumber = contignumber;
  }

  return 0;
}

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
