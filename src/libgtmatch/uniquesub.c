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

#include <inttypes.h>
#include <string.h>
#include "libgtcore/error.h"
#include "libgtcore/seqiterator.h"
#include "defined-types.h"
#include "alphadef.h"
#include "spacedef.h"
#include "optionargmode.h"
#include "format64.h"
#include "uniquesub.h"
#include "encseq-def.h"

typedef struct
{
  bool showsequence,
       showquerypos,
       showsubjectpos;
  Definedunsignedlong minlength,
                      maxlength;
} Rangespecinfo;

typedef int (*Preprocessuniquelength)(uint64_t,
                                      const char *,
                                      void *,
                                      Error *err);
typedef int (*Processuniquelength)(const Alphabet *,
                                   const Uchar *,
                                   unsigned long,
                                   unsigned long,
                                   Seqpos,
                                   void *,
                                   Error *err);
typedef int (*Postprocessuniquelength)(const Alphabet *,
                                       uint64_t,
                                       const char *,
                                       const Uchar *,
                                       unsigned long,
                                       void *,
                                       Error *err);

typedef struct
{
  const void *genericindex;
  const Alphabet *alphabet;
  Uniqueforwardfunction uniqueforward;
  Preprocessuniquelength preprocessuniquelength;
  Processuniquelength processuniquelength;
  Postprocessuniquelength postprocessuniquelength;
  void *processinfo;
  const Encodedsequence *encseq;
} Substringinfo;

static void checkifsequenceisthere(const Encodedsequence *encseq,
                                   Seqpos witnessposition,
                                   unsigned long uniquelength,
                                   const Uchar *qptr)
{
  unsigned long i;

  for (i=0; i<uniquelength; i++)
  {
    if (qptr[i] != getencodedchar(encseq,witnessposition+i,Forwardmode))
    {
      fprintf(stderr,"at witnesspos " FormatSeqpos " query[%lu] = %u != %u = "
                     " subject[%lu]\n",
                     PRINTSeqposcast(witnessposition), 
                     i,
                     qptr[i],
                     getencodedchar(encseq,witnessposition+i,Forwardmode),
                     witnessposition+i);
      exit(EXIT_FAILURE); /* Program error */
    }
  }
}

static unsigned long checksperformed = 0;

static int uniqueposinsinglesequence(Substringinfo *substringinfo,
                                     uint64_t unitnum,
                                     const Uchar *query,
                                     unsigned long querylen,
                                     const char *desc,
                                     Error *err)
{
  const Uchar *qptr;
  unsigned long uniquelength, remaining;
  Seqpos witnessposition, *wptr;

  error_check(err);
  if (substringinfo->preprocessuniquelength != NULL &&
      substringinfo->preprocessuniquelength(unitnum,
                                            desc,
                                            substringinfo->processinfo,
                                            err) != 0)
  {
    return -1;
  }
  if (((Rangespecinfo *) substringinfo->processinfo)->showsubjectpos ||
      substringinfo->encseq != NULL)
  {
    wptr = &witnessposition;
  } else
  {
    wptr = NULL;
  }
  for (qptr = query, remaining = querylen; remaining > 0; qptr++, remaining--)
  {
    uniquelength = substringinfo->uniqueforward(substringinfo->genericindex,
                                                wptr,
                                                qptr,
                                                query+querylen);
    if (uniquelength > 0)
    {
      if (wptr != NULL && substringinfo->encseq != NULL)
      {
        checksperformed++;
        checkifsequenceisthere(substringinfo->encseq,
                               witnessposition,
                               uniquelength,
                               qptr);
      }
      if (substringinfo->processuniquelength(substringinfo->alphabet,
                                             query,
                                             uniquelength,
                                             (unsigned long) (qptr-query),
                                             wptr == NULL
                                              ? (Seqpos) 0
                                              : witnessposition,
                                             substringinfo->processinfo,
                                             err) != 0)
      {
        return -2;
      }
    }
  }
  if (substringinfo->postprocessuniquelength != NULL &&
      substringinfo->postprocessuniquelength(substringinfo->alphabet,
                                             unitnum,
                                             desc,
                                             query,
                                             querylen,
                                             substringinfo->processinfo,
                                             err) != 0)
  {
    return -3;
  }
  return 0;
}

static int showunitnum(uint64_t unitnum,
                       const char *desc,
                       /*@unused@*/ void *info,
                       /*@unused@*/ Error *err)
{
  printf("unit " Formatuint64_t, PRINTuint64_tcast(unitnum));
  if (desc != NULL && desc[0] != '\0')
  {
    printf(" (%s)",desc);
  }
  printf("\n");
  return 0;
}

static int showifinlengthrange(const Alphabet *alphabet,
                               const Uchar *start,
                               unsigned long uniquelength,
                               unsigned long querystart,
                               Seqpos subjectpos,
                               void *info,
                                /*@unused@*/ Error *err)
{
  Rangespecinfo *rangespecinfo = (Rangespecinfo *) info;

  if ((!rangespecinfo->minlength.defined ||
      uniquelength >= rangespecinfo->minlength.valueunsignedlong) &&
     (!rangespecinfo->maxlength.defined ||
      uniquelength <= rangespecinfo->maxlength.valueunsignedlong))
  {
    if (rangespecinfo->showquerypos)
    {
      printf("%lu ",querystart);
    }
    printf("%lu",uniquelength);
    if (rangespecinfo->showsubjectpos)
    {
      printf(" " FormatSeqpos,PRINTSeqposcast(subjectpos));
    }
    if (rangespecinfo->showsequence)
    {
      (void) putchar(' ');
      showsymbolstring(alphabet,
                       start + querystart,
                       uniquelength);
    }
    (void) putchar('\n');
  }
  return 0;
}

int findsubqueryuniqueforward(const Encodedsequence *encseq,
                              const void *genericindex,
                              Uniqueforwardfunction uniqueforward,
                              const Alphabet *alphabet,
                              const StrArray *queryfilenames,
                              Definedunsignedlong minlength,
                              Definedunsignedlong maxlength,
                              bool showsequence,
                              bool showquerypos,
                              bool showsubjectpos,
                              Error *err)
{
  Substringinfo substringinfo;
  Rangespecinfo rangespecinfo;
  bool haserr = false;
  SeqIterator *seqit;
  const Uchar *query;
  unsigned long querylen;
  char *desc = NULL;
  int retval;
  uint64_t unitnum;

  error_check(err);
  substringinfo.genericindex = genericindex;
  rangespecinfo.minlength = minlength;
  rangespecinfo.maxlength = maxlength;
  rangespecinfo.showsequence = showsequence;
  rangespecinfo.showquerypos = showquerypos;
  rangespecinfo.showsubjectpos = showsubjectpos;
  substringinfo.preprocessuniquelength = showunitnum;
  substringinfo.processuniquelength = showifinlengthrange;
  substringinfo.postprocessuniquelength = NULL;
  substringinfo.alphabet = alphabet;
  substringinfo.processinfo = &rangespecinfo;
  substringinfo.uniqueforward = uniqueforward;
  substringinfo.encseq = encseq;
  seqit = seqiterator_new(queryfilenames,getsymbolmapAlphabet(alphabet),true);
  for (unitnum = 0; /* Nothing */; unitnum++)
  {
    retval = seqiterator_next(seqit,
                              &query,
                              &querylen,
                              &desc,
                              err);
    if (retval < 0)
    {
      haserr = true;
      break;
    }
    if (retval == 0)
    {
      break;
    }
    if (uniqueposinsinglesequence(&substringinfo,
                                  unitnum,
                                  query,
                                  querylen,
                                  desc,
                                  err) != 0)
    {
      haserr = true;
      break;
    }
    FREESPACE(desc);
  }
  seqiterator_delete(seqit);
  /* printf("# %lu checks performed\n",checksperformed); */
  return haserr ? -1 : 0;
}
