/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdio.h>
#include <stdlib.h>
#include "libgtcore/env.h"
#include "libgtcore/minmax.h"
#include "types.h"
#include "chardef.h"
#include "intbits-tab.h"
#include "encseq-def.h"

#include "sfx-cmpsuf.pr"

static void showlocalsuffix(FILE *fpout,
                            const Encodedsequence *encseq,
                            const Uchar *characters,
                            Seqpos start,
                            Seqpos depth)
{
  Seqpos i, end, totallength = getencseqtotallength(encseq);
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
  showlocalsuffix(stderr,encseq,characters,*ptr1,depth);
  fprintf(stderr,"\",\"");
  showlocalsuffix(stderr,encseq,characters,*ptr2,depth);
  fprintf(stderr,"\"=" FormatSeqpos ")=%d with maxlcp " FormatSeqpos "\n",
              PRINTSeqposcast(*ptr2),
              cmp,
              PRINTSeqposcast(maxlcp));
}

void checkifprefixesareidentical(const Encodedsequence *encseq,
                                 const Uchar *characters,
                                 const Seqpos *suftab,
                                 uint32_t prefixlength,
                                 Seqpos depth,
                                 Seqpos left,
                                 Seqpos right)
{
  const Seqpos *ptr;
  Seqpos maxlcp;
  int cmp;

  for (ptr = suftab + left; ptr < suftab + right; ptr++)
  {
    cmp = comparetwosuffixes(encseq,
                             &maxlcp,
                             false,
                             true,
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
                            ptr,ptr+1,cmp,maxlcp);
      exit(EXIT_FAILURE);
    }
  }
}

void showentiresuftab(const Encodedsequence *encseq,
                      const Uchar *characters,
                      const Seqpos *suftab,
                      Seqpos depth)
{
  const Seqpos *ptr;
  Seqpos totallength = getencseqtotallength(encseq);

  for (ptr = suftab; ptr <= suftab + totallength; ptr++)
  {
    printf("suftab[" FormatSeqpos "]=" FormatSeqpos " ",
            PRINTSeqposcast((Seqpos) (ptr-suftab)),
            PRINTSeqposcast(*ptr));
    showlocalsuffix(stdout,encseq,characters,*ptr,depth);
    printf("\n");
  }
}

void checkentiresuftab(const Encodedsequence *encseq,
                       const Uchar *characters,
                       const Seqpos *suftab,
                       bool specialsareequal,
                       bool specialsareequalatdepth0,
                       Seqpos depth,
                       Env *env)
{
  const Seqpos *ptr;
  Bitstring *startposoccurs;
  Seqpos maxlcp, countbitsset = 0, totallength = getencseqtotallength(encseq);
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
    cmp = comparetwosuffixes(encseq,
                             &maxlcp,
                             specialsareequal,
                             specialsareequalatdepth0,
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
