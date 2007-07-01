/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdio.h>
#include <stdlib.h>
#include "libgtcore/env.h"
#include "types.h"
#include "chardef.h"
#include "minmax.h"
#include "intbits-tab.h"
#include "encseq-def.h"

static int comparetwoprefixes(const Encodedsequence *encseq,
                              Seqpos *maxlcp,
                              bool specialsareequal,
                              bool specialsareequalatdepth0,
                              Seqpos totallength,
                              Seqpos depth,
                              Seqpos start1,
                              Seqpos start2)
{
  Uchar cc1, cc2;
  Seqpos pos1, pos2, end1, end2;

  end1 = end2 = totallength;
  if (depth > 0)
  {
    if (end1 > start1 + depth)
    {
      end1 = start1 + depth;
    }
    if (end2 > start2 + depth)
    {
      end2 = start2 + depth;
    }
  }
  for (pos1=start1, pos2=start2; pos1 < end1 && pos2 < end2; pos1++, pos2++)
  {
    cc1 = getencodedchar(encseq,pos1);
    cc2 = getencodedchar(encseq,pos2);
    if (ISSPECIAL(cc1))
    {
      if (ISSPECIAL(cc2))
      {
        if (specialsareequal || (pos1 == start1 && specialsareequalatdepth0))
        {
          *maxlcp = pos1 - start1 + 1;
          return 0;
        }
        if (pos1 < pos2)
        {
          *maxlcp = pos1  - start1;
          return -1; /* a < b */
        }
        if (pos1 > pos2)
        {
          *maxlcp = pos1 - start1;
          return 1; /* a > b */
        }
        *maxlcp = pos1 - start1 + 1;
        return 0; /* a = b */
      }
      *maxlcp = pos1 - start1;
      return 1; /* a > b */
    } else
    {
      if (ISSPECIAL(cc2))
      {
        *maxlcp = pos1 - start1;
        return -1; /* a < b */
      }
      if (cc1 < cc2)
      {
        *maxlcp = pos1 - start1;
        return -1;
      }
      if (cc1 > cc2)
      {
        *maxlcp = pos1 - start1;
        return 1;
      }
    }
  }
  *maxlcp = pos1 - start1;
  return 0;
}

static void showlocalsuffix(FILE *fpout,
                            const Encodedsequence *encseq,
                            const Uchar *characters,
                            Seqpos start,
                            Seqpos depth,
                            Seqpos totallength)
{
  Seqpos i, end;
  Uchar cc;
  const Seqpos maxshow = (Seqpos) 30;

  if (depth == 0)
  {
    end = MIN(maxshow,totallength);
  } else
  {
    end = MIN(maxshow,MIN(totallength,start+depth-1));
  }
  for (i = start; i <= end; i++)
  {
    if (i == totallength)
    {
      (void) putc('~',fpout);
      break;
    }
    cc = getencodedchar(encseq,i);
    if (ISSPECIAL(cc))
    {
      (void) putc('~',fpout);
      break;
    }
    (void) putc((Fputcfirstargtype) characters[(int) cc],fpout);
  }
}

static void showcomparisonfailure(const char *where,
                                  const Encodedsequence *encseq,
                                  const Uchar *characters,
                                  const Seqpos *suftab,
                                  Seqpos depth,
                                  Seqpos totallength,
                                  const Seqpos *ptr1,
                                  const Seqpos *ptr2,
                                  int cmp,
                                  Seqpos maxlcp)
{
  fprintf(stderr,"ERROR: %s(" FormatSeqpos " vs " FormatSeqpos 
                 " " FormatSeqpos "=\"",
                       where,
                       PRINTSeqposcast((Seqpos) (ptr1 - suftab)),
                       PRINTSeqposcast((Seqpos) (ptr2 - suftab)),
                       PRINTSeqposcast(*ptr1));
  showlocalsuffix(stderr,encseq,characters,*ptr1,depth,totallength);
  fprintf(stderr,"\",\"");
  showlocalsuffix(stderr,encseq,characters,*ptr2,depth,totallength);
  fprintf(stderr,"\"=" FormatSeqpos ")=%d with maxlcp " FormatSeqpos "\n",
              PRINTSeqposcast(*ptr2),
              cmp,
              PRINTSeqposcast(maxlcp));
}

void checkifprefixesareidentical(const Encodedsequence *encseq,
                                 const Uchar *characters,
                                 const Seqpos *suftab,
                                 uint32_t prefixlength,
                                 Seqpos totallength,
                                 Seqpos depth,
                                 Seqpos left,
                                 Seqpos right)
{
  const Seqpos *ptr;
  Seqpos maxlcp;
  int cmp;

  for (ptr = suftab + left; ptr < suftab + right; ptr++)
  {
    cmp = comparetwoprefixes(encseq,
                             &maxlcp,
                             false,
                             true,
                             totallength,
                             depth,
                             *ptr,
                             *(ptr+1));
    if (cmp != 0 || maxlcp != (Seqpos) prefixlength)
    {
      showcomparisonfailure("checkifprefixesareidentical",
                            encseq,
                            characters,
                            suftab,
                            depth,
                            totallength,
                            ptr,ptr+1,cmp,maxlcp);
      exit(EXIT_FAILURE);
    }
  }
}

void showentiresuftab(const Encodedsequence *encseq,
                      const Uchar *characters,
                      const Seqpos *suftab,
                      Seqpos depth,
                      Seqpos totallength)
{
  const Seqpos *ptr;

  for (ptr = suftab; ptr <= suftab + totallength; ptr++)
  {
    printf("suftab[" FormatSeqpos "]=" FormatSeqpos " ",
            PRINTSeqposcast((Seqpos) (ptr-suftab)),
            PRINTSeqposcast(*ptr));
    showlocalsuffix(stdout,encseq,characters,*ptr,depth,totallength);
    printf("\n");
  }
}

void checkentiresuftab(const Encodedsequence *encseq,
                       const Uchar *characters,
                       const Seqpos *suftab,
                       bool specialsareequal,
                       bool specialsareequalatdepth0,
                       Seqpos totallength,
                       Seqpos depth,
                       Env *env)
{
  const Seqpos *ptr;
  Bitstring *startposoccurs;
  Seqpos maxlcp, countbitsset = 0;
  int cmp;

  assert(!specialsareequal || specialsareequalatdepth0);
  INITBITTAB(startposoccurs,totallength+1);
  for (ptr = suftab; ptr <= suftab + totallength; ptr++)
  {
    if (ISIBITSET(startposoccurs,*ptr))
    {
      fprintf(stderr,"ERROR: suffix with startpos " FormatSeqpos 
                     " already occurs\n",
                      PRINTSeqposcast(*ptr));
      exit(EXIT_FAILURE);
    }
    SETIBIT(startposoccurs,*ptr);
    countbitsset++;
  }
  if (countbitsset != totallength+1)
  {
    fprintf(stderr,"ERROR: not all bits are set\n");
    exit(EXIT_FAILURE);
  }
  FREESPACE(startposoccurs);
  for (ptr = suftab + 1; ptr <= suftab + totallength; ptr++)
  {
    cmp = comparetwoprefixes(encseq,
                             &maxlcp,
                             specialsareequal,
                             specialsareequalatdepth0,
                             totallength,
                             depth,
                             *(ptr-1),
                             *ptr);
    if (cmp > 0)
    {
      showcomparisonfailure("checkentiresuftab",
                            encseq,
                            characters,
                            suftab,
                            depth,
                            totallength,
                            ptr-1,
                            ptr,
                            cmp,
                            maxlcp);
      exit(EXIT_FAILURE);
    }
  }
  /*
  printf("# checkentiresuftab with mode 'specials are %s'\n",
               specialsareequal ? "equal" : "different");
  */
}
