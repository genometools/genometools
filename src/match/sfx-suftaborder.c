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
#include "core/error_api.h"
#include "core/minmax.h"
#include "core/intbits.h"
#include "core/encseq.h"
#include "esa-seqread.h"
#include "sfx-suftaborder.h"
#include "suffixptr.h"

static void showlocalsuffix(FILE *fpout,
                            const GtEncseq *encseq,
                            GtReadmode readmode,
                            unsigned long start,
                            unsigned long depth)
{
  unsigned long i, end, totallength;
  GtUchar cc;
  const unsigned long maxshow = (unsigned long) 30;
  const GtUchar *characters;

  totallength = gt_encseq_total_length(encseq);
  characters = gt_alphabet_characters(gt_encseq_alphabet(encseq));
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
    cc = gt_encseq_get_encoded_char(encseq,i,readmode);
    if (ISSPECIAL(cc))
    {
      (void) putc('~',fpout);
      break;
    }
    (void) putc((int) characters[(int) cc],fpout);
  }
}

static void showcomparisonfailure(const char *filename,
                                  int line,
                                  const char *where,
                                  const GtEncseq *encseq,
                                  GtReadmode readmode,
                                  const Suffixptr *suftab,
                                  unsigned long depth,
                                  const Suffixptr *ptr1,
                                  const Suffixptr *ptr2,
                                  int cmp,
                                  unsigned long maxlcp)
{
  fprintf(stderr,"ERROR: file \"%s\", line %d: ",filename,line);
  fprintf(stderr,"%s(%lu vs %lu"
                 " %lu=\"",
                       where,
                       (unsigned long) (ptr1 - suftab),
                       (unsigned long) (ptr2 - suftab),
                       SUFFIXPTRDEREF(ptr1));
  showlocalsuffix(stderr,encseq,readmode,SUFFIXPTRDEREF(ptr1),depth);
  fprintf(stderr,"\",\"");
  showlocalsuffix(stderr,encseq,readmode,SUFFIXPTRDEREF(ptr2),depth);
  fprintf(stderr,"\"=%lu)=%d with maxlcp %lu\n",SUFFIXPTRDEREF(ptr2),
                                                cmp,
                                                maxlcp);
}

void gt_checkifprefixesareidentical(const char *filename,
                                 int line,
                                 const GtEncseq *encseq,
                                 GtReadmode readmode,
                                 const Suffixptr *suftab,
                                 unsigned int prefixlength,
                                 unsigned long depth,
                                 unsigned long left,
                                 unsigned long right)
{
  const Suffixptr *ptr;
  unsigned long maxlcp;
  int cmp;
  GtEncseqReader *esr1, *esr2;
  bool haserr = false;

  esr1 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  esr2 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  for (ptr = suftab + left; ptr < suftab + right; ptr++)
  {
    cmp = gt_encseq_comparetwosuffixes(encseq,
                                       readmode,
                                       &maxlcp,
                                       false,
                                       true,
                                       depth,
                                       SUFFIXPTRDEREF(ptr),
                                       SUFFIXPTRDEREF(ptr+1),
                                       esr1,
                                       esr2);
    if (cmp != 0 || maxlcp != (unsigned long) prefixlength)
    {
      showcomparisonfailure(filename,
                            line,
                            "checkifprefixesareidentical",
                            encseq,
                            readmode,
                            suftab,
                            depth,
                            ptr,ptr+1,cmp,maxlcp);
      haserr = true;
      break;
    }
  }
  gt_encseq_reader_delete(esr1);
  gt_encseq_reader_delete(esr2);
  if (haserr)
  {
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

void gt_showentiresuftab(const GtEncseq *encseq,
                         GtReadmode readmode,
                         const Suffixptr *suftab,
                         unsigned long depth)
{
  const Suffixptr *ptr;
  unsigned long totallength = gt_encseq_total_length(encseq);

  for (ptr = suftab; ptr <= suftab + totallength; ptr++)
  {
    printf("suftab[%lu]=%lu ",
            (unsigned long) (ptr-suftab),
            SUFFIXPTRDEREF(ptr));
    showlocalsuffix(stdout,encseq,readmode,SUFFIXPTRDEREF(ptr),depth);
    printf("\n");
  }
}

void gt_checksortedsuffixes(const char *filename,
                            int line,
                            const GtEncseq *encseq,
                            GtReadmode readmode,
                            const Suffixptr *suftab,
                            unsigned long numberofsuffixes,
                            bool specialsareequal,
                            bool specialsareequalatdepth0,
                            unsigned long depth)
{
  const Suffixptr *ptr;
  unsigned long maxlcp, totallength = gt_encseq_total_length(encseq);
  GtEncseqReader *esr1, *esr2;
  int cmp;

  gt_assert(!specialsareequal || specialsareequalatdepth0);
  esr1 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  esr2 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  gt_assert(numberofsuffixes > 0);
  gt_assert(SUFFIXPTRDEREF(suftab) < totallength);
  for (ptr = suftab + 1; ptr < suftab + numberofsuffixes; ptr++)
  {
    if (ptr < suftab + numberofsuffixes - 1)
    {
      gt_assert(SUFFIXPTRDEREF(ptr) < totallength);
      cmp = gt_encseq_comparetwosuffixes(encseq,
                                         readmode,
                                         &maxlcp,
                                         specialsareequal,
                                         specialsareequalatdepth0,
                                         depth,
                                         SUFFIXPTRDEREF(ptr-1),
                                         SUFFIXPTRDEREF(ptr),
                                         esr1,
                                         esr2);
      if (cmp > 0)
      {
        showcomparisonfailure(filename,
                              line,
                              "checksortedsuffixes",
                              encseq,
                              readmode,
                              suftab,
                              depth,
                              ptr-1,
                              ptr,
                              cmp,
                              maxlcp);
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
    } else
    {
      if (numberofsuffixes == totallength+1)
      {
        gt_assert(SUFFIXPTRDEREF(ptr) == totallength);
      }
    }
  }
  gt_encseq_reader_delete(esr1);
  gt_encseq_reader_delete(esr2);
}

void gt_checkentiresuftab(const char *filename,
                          int line,
                          const GtEncseq *encseq,
                          GtReadmode readmode,
                          const Suffixptr *suftab,
                          unsigned long numberofsuffixes,
                          Sequentialsuffixarrayreader *ssar,
                          bool specialsareequal,
                          bool specialsareequalatdepth0,
                          unsigned long depth,
                          GtError *err)
{
  const Suffixptr *ptr;
  unsigned long maxlcp,
                currentlcp = 0,
                totallength = gt_encseq_total_length(encseq);
  int cmp;
  GtEncseqReader *esr1, *esr2;
  bool haserr = false;

#ifdef INLINEDSequentialsuffixarrayreader
  GtUchar tmpsmalllcpvalue;
#else
  int retval;
#endif

  gt_error_check(err);
  gt_assert(!specialsareequal || specialsareequalatdepth0);
  if (numberofsuffixes == totallength+1)
  {
    GtBitsequence *startposoccurs;
    unsigned long countbitsset = 0;

    GT_INITBITTAB(startposoccurs,totallength+1);
    for (ptr = suftab; ptr <= suftab + totallength; ptr++)
    {
      if (GT_ISIBITSET(startposoccurs,SUFFIXPTRDEREF(ptr)))
      {
        fprintf(stderr,"ERROR: suffix with startpos %lu"
                       " already occurs\n",
                        SUFFIXPTRDEREF(ptr));
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
      GT_SETIBIT(startposoccurs,SUFFIXPTRDEREF(ptr));
      countbitsset++;
    }
    if (countbitsset != totallength+1)
    {
      fprintf(stderr,"ERROR: not all bits are set\n");
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    gt_free(startposoccurs);
  }
  esr1 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  esr2 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  gt_assert(numberofsuffixes > 0);
  gt_assert(SUFFIXPTRDEREF(suftab) < totallength);
  for (ptr = suftab + 1; !haserr && ptr < suftab + numberofsuffixes; ptr++)
  {
    if (ptr < suftab + numberofsuffixes - 1)
    {
      gt_assert(SUFFIXPTRDEREF(ptr) < totallength);
      cmp = gt_encseq_comparetwosuffixes(encseq,
                                         readmode,
                                         &maxlcp,
                                         specialsareequal,
                                         specialsareequalatdepth0,
                                         depth,
                                         SUFFIXPTRDEREF(ptr-1),
                                         SUFFIXPTRDEREF(ptr),
                                         esr1,
                                         esr2);
      if (cmp > 0)
      {
        showcomparisonfailure(filename,
                              line,
                              "checkentiresuftab",
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
      if (numberofsuffixes == totallength+1)
      {
        gt_assert(SUFFIXPTRDEREF(ptr) == totallength);
      }
    }
    if (ssar != NULL)
    {
#ifdef INLINEDSequentialsuffixarrayreader
      NEXTSEQUENTIALLCPTABVALUE(currentlcp,ssar);
#else
      retval = gt_nextSequentiallcpvalue(&currentlcp,ssar,err);
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
        fprintf(stderr,"%lu: startpos=%lu, firstchar=%u, "
                "startpos=%lu,firstchar=%u",
                (unsigned long) (ptr - suftab),
                SUFFIXPTRDEREF(ptr-1),
                (unsigned int) gt_encseq_get_encoded_char(encseq,
                                                          SUFFIXPTRDEREF(ptr-1),
                                                          readmode),
                SUFFIXPTRDEREF(ptr),
                (SUFFIXPTRDEREF(ptr) < totallength)
                   ? (unsigned int) gt_encseq_get_encoded_char(encseq,
                                                            SUFFIXPTRDEREF(ptr),
                                                            readmode)
                   : SEPARATOR);
        fprintf(stderr,", maxlcp(bruteforce) = %lu != %lu(fast)\n",
                          maxlcp, currentlcp);
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
    }
  }
  gt_encseq_reader_delete(esr1);
  gt_encseq_reader_delete(esr2);
  if (haserr)
  {
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  /*
  printf("# gt_checkentiresuftab with mode 'specials are %s'\n",
               specialsareequal ? "equal" : "different");
  */
}
