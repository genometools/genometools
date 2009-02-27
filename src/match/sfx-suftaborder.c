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
#include "core/chardef.h"
#include "core/error.h"
#include "core/minmax.h"
#include "intbits-tab.h"
#include "encseq-def.h"
#include "esa-seqread.h"

#include "sfx-cmpsuf.pr"

static void showlocalsuffix(FILE *fpout,
                            const Encodedsequence *encseq,
                            Readmode readmode,
                            Seqpos start,
                            Seqpos depth)
{
  Seqpos i, end, totallength;
  Uchar cc;
  const Seqpos maxshow = (Seqpos) 30;
  const Uchar *characters;

  totallength = getencseqtotallength(encseq);
  characters = getencseqAlphabetcharacters(encseq);
  if (depth == 0)
  {
    end = MIN(start + maxshow,totallength);
  } else
  {
    end = MIN(start + maxshow,MIN(totallength,start+depth));
  }
  for (i = start; i <= end; i++)
  {
    if (i == totallength)
    {
      (void) putc('~',fpout);
      break;
    }
    cc = getencodedchar(encseq,i,readmode); /* for testing */
    if (ISSPECIAL(cc))
    {
      (void) putc('~',fpout);
      break;
    }
    (void) putc((int) characters[(int) cc],fpout);
  }
}

static void showcomparisonfailure(const char *where,
                                  const Encodedsequence *encseq,
                                  Readmode readmode,
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
  showlocalsuffix(stderr,encseq,readmode,*ptr1,depth);
  fprintf(stderr,"\",\"");
  showlocalsuffix(stderr,encseq,readmode,*ptr2,depth);
  fprintf(stderr,"\"=" FormatSeqpos ")=%d with maxlcp " FormatSeqpos "\n",
              PRINTSeqposcast(*ptr2),
              cmp,
              PRINTSeqposcast(maxlcp));
}

void checkifprefixesareidentical(const Encodedsequence *encseq,
                                 Readmode readmode,
                                 const Seqpos *suftab,
                                 unsigned int prefixlength,
                                 Seqpos depth,
                                 Seqpos left,
                                 Seqpos right)
{
  const Seqpos *ptr;
  Seqpos maxlcp;
  int cmp;
  Encodedsequencescanstate *esr1, *esr2;
  bool haserr = false;

  esr1 = newEncodedsequencescanstate();
  esr2 = newEncodedsequencescanstate();
  for (ptr = suftab + left; ptr < suftab + right; ptr++)
  {
    cmp = comparetwosuffixes(encseq,
                             readmode,
                             &maxlcp,
                             false,
                             true,
                             depth,
                             *ptr,
                             *(ptr+1),
                             esr1,
                             esr2);
    if (cmp != 0 || maxlcp != (Seqpos) prefixlength)
    {
      showcomparisonfailure("checkifprefixesareidentical",
                            encseq,
                            readmode,
                            suftab,
                            depth,
                            ptr,ptr+1,cmp,maxlcp);
      haserr = true;
      break;
    }
  }
  freeEncodedsequencescanstate(&esr1);
  freeEncodedsequencescanstate(&esr2);
  if (haserr)
  {
    exit(EXIT_FAILURE); /* programming error */
  }
}

void showentiresuftab(const Encodedsequence *encseq,
                      Readmode readmode,
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
    showlocalsuffix(stdout,encseq,readmode,*ptr,depth);
    printf("\n");
  }
}

void checkentiresuftab(const Encodedsequence *encseq,
                       Readmode readmode,
                       const Seqpos *suftab,
                       Sequentialsuffixarrayreader *ssar,
                       bool specialsareequal,
                       bool specialsareequalatdepth0,
                       Seqpos depth,
                       GtError *err)
{
  const Seqpos *ptr;
  Bitstring *startposoccurs;
  Seqpos maxlcp, countbitsset = 0, currentlcp = 0,
         totallength = getencseqtotallength(encseq);
  int cmp;
  Encodedsequencescanstate *esr1, *esr2;
  bool haserr = false;

#ifdef INLINEDSequentialsuffixarrayreader
  Uchar tmpsmalllcpvalue;
#else
  int retval;
#endif

  gt_error_check(err);
  gt_assert(!specialsareequal || specialsareequalatdepth0);
  INITBITTAB(startposoccurs,totallength+1);
  for (ptr = suftab; ptr <= suftab + totallength; ptr++)
  {
    if (ISIBITSET(startposoccurs,*ptr))
    {
      fprintf(stderr,"ERROR: suffix with startpos " FormatSeqpos
                     " already occurs\n",
                      PRINTSeqposcast(*ptr));
      exit(EXIT_FAILURE); /* programming error */
    }
    SETIBIT(startposoccurs,*ptr);
    countbitsset++;
  }
  if (countbitsset != totallength+1)
  {
    fprintf(stderr,"ERROR: not all bits are set\n");
    exit(EXIT_FAILURE); /* programming error */
  }
  gt_free(startposoccurs);
  esr1 = newEncodedsequencescanstate();
  esr2 = newEncodedsequencescanstate();
  gt_assert(*suftab < totallength);
  for (ptr = suftab + 1; !haserr && ptr <= suftab + totallength; ptr++)
  {
    if (ptr < suftab + totallength)
    {
      gt_assert(*ptr < totallength);
      cmp = comparetwosuffixes(encseq,
                               readmode,
                               &maxlcp,
                               specialsareequal,
                               specialsareequalatdepth0,
                               depth,
                               *(ptr-1),
                               *ptr,
                               esr1,
                               esr2);
      if (cmp > 0)
      {
        showcomparisonfailure("checkentiresuftab",
                              encseq,
                              readmode,
                              suftab,
                              depth,
                              ptr-1,
                              ptr,
                              cmp,
                              maxlcp);
        haserr = true;
        break;
      }
    } else
    {
      maxlcp = 0;
      gt_assert(*ptr == totallength);
    }
    if (ssar != NULL)
    {
#ifdef INLINEDSequentialsuffixarrayreader
      NEXTSEQUENTIALLCPTABVALUE(currentlcp,ssar);
#else
      retval = nextSequentiallcpvalue(&currentlcp,ssar,err);
      if (retval < 0)
      {
        haserr = true;
        break;
      }
      if (retval == 0)
      {
        break;
      }
#endif
      if (maxlcp != currentlcp)
      {
        fprintf(stderr,"%lu: startpos=" FormatSeqpos ", firstchar=%u, "
                "startpos=" FormatSeqpos ",firstchar=%u",
                (unsigned long) (ptr - suftab),
                PRINTSeqposcast(*(ptr-1)),
                (unsigned int) getencodedchar(encseq,*(ptr-1),readmode),
                PRINTSeqposcast(*ptr),
                (*ptr < totallength)
                ? (unsigned int) getencodedchar(encseq,*ptr,readmode)
                : SEPARATOR);
        fprintf(stderr,", maxlcp(bruteforce) = " FormatSeqpos " != "
                          FormatSeqpos "(fast)\n",
                    PRINTSeqposcast(maxlcp),
                    PRINTSeqposcast(currentlcp));
        /* exit(EXIT_FAILURE); programming error */
      }
    }
  }
  freeEncodedsequencescanstate(&esr1);
  freeEncodedsequencescanstate(&esr2);
  if (haserr)
  {
    exit(EXIT_FAILURE); /* programming error */
  }
  /*
  printf("# checkentiresuftab with mode 'specials are %s'\n",
               specialsareequal ? "equal" : "different");
  */
}
