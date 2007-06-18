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
#include "intbits.h"
#include "encseq-def.h"

static int comparetwoprefixes(const Encodedsequence *encseq,
                              Uint *maxlcp,
                              bool specialsareequal,
                              bool specialsareequalatdepth0,
                              Uint totallength,
                              Uint depth,
                              Uint start1,
                              Uint start2)
{
  Uchar cc1, cc2;
  Uint pos1, pos2, end1, end2;

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
                            Uint start,
                            Uint depth,
                            Uint totallength)
{
  Uint i, end;
  Uchar cc;
  const Uint maxshow = UintConst(30);

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
                                  const Uint *suftab,
                                  Uint depth,
                                  Uint totallength,
                                  const Uint *ptr1,
                                  const Uint *ptr2,
                                  int cmp,
                                  Uint maxlcp)
{
  fprintf(stderr,"ERROR: %s(%lu vs %lu: %lu=\"",
                       where,
                       (Showuint) (ptr1 - suftab),
                       (Showuint) (ptr2 - suftab),
                       (Showuint) *ptr1);
  showlocalsuffix(stderr,encseq,characters,*ptr1,depth,totallength);
  fprintf(stderr,"\",\"");
  showlocalsuffix(stderr,encseq,characters,*ptr2,depth,totallength);
  fprintf(stderr,"\"=%lu)=%d with maxlcp %lu\n",
              (Showuint) *ptr2,cmp,(Showuint) maxlcp);
}

void checkifprefixesareidentical(const Encodedsequence *encseq,
                                 const Uchar *characters,
                                 const Uint *suftab,
                                 unsigned int prefixlength,
                                 Uint64 totallength,
                                 Uint depth,
                                 Uint left,
                                 Uint right)
{
  const Uint *ptr;
  Uint maxlcp;
  int cmp;

  CHECKIFFITS32BITS(totallength);
  for (ptr = suftab + left; ptr < suftab + right; ptr++)
  {
    cmp = comparetwoprefixes(encseq,
                             &maxlcp,
                             false,
                             true,
                             (Uint) totallength,
                             depth,
                             *ptr,
                             *(ptr+1));
    if (cmp != 0 || prefixlength != maxlcp)
    {
      showcomparisonfailure("checkifprefixesareidentical",
                            encseq,
                            characters,
                            suftab,
                            depth,
                            (Uint) totallength,
                            ptr,ptr+1,cmp,maxlcp);
      exit(EXIT_FAILURE);
    }
  }
}

void showentiresuftab(const Encodedsequence *encseq,
                      const Uchar *characters,
                      const Uint *suftab,
                      Uint depth,
                      Uint totallength)
{
  const Uint *ptr;

  for (ptr = suftab; ptr <= suftab + totallength; ptr++)
  {
    printf("suftab[%lu]=%lu ",(Showuint) (ptr-suftab),(Showuint) *ptr);
    showlocalsuffix(stdout,encseq,characters,*ptr,depth,totallength);
    printf("\n");
  }
}

void checkentiresuftab(const Encodedsequence *encseq,
                       const Uchar *characters,
                       const Uint *suftab,
                       bool specialsareequal,
                       bool specialsareequalatdepth0,
                       Uint64 totallength,
                       Uint depth,
                       Env *env)
{
  const Uint *ptr;
  Uint *startposoccurs, maxlcp, countbitsset = 0;
  int cmp;

  assert(!specialsareequal || specialsareequalatdepth0);
  CHECKIFFITS32BITS(totallength+1);
  INITBITTAB(startposoccurs,(Uint) (totallength+1));
  for (ptr = suftab; ptr <= suftab + totallength; ptr++)
  {
    if (ISIBITSET(startposoccurs,*ptr))
    {
      fprintf(stderr,"ERROR: suffix with startpos %lu already occurs\n",
                      (Showuint) *ptr);
      exit(EXIT_FAILURE);
    }
    SETIBIT(startposoccurs,*ptr);
    countbitsset++;
  }
  if (countbitsset != (Uint) (totallength+1))
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
                             (Uint) totallength,
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
                            (Uint) totallength,
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
