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

typedef struct
{
  bool showsequence,
       showquerypos,
       showrefpos;
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
} Substringinfo;

static int uniqueposinsinglesequence(Substringinfo *substringinfo,
                                     uint64_t unitnum,
                                     const Uchar *query,
                                     unsigned long querylen,
                                     const char *desc,
                                     Error *err)
{
  const Uchar *qptr;
  unsigned long uniquelength, remaining;
  Seqpos witnessposition;

  error_check(err);
  if (substringinfo->preprocessuniquelength != NULL &&
      substringinfo->preprocessuniquelength(unitnum,
                                            desc,
                                            substringinfo->processinfo,
                                            err) != 0)
  {
    return -1;
  }
  for (qptr = query, remaining = querylen; remaining > 0; qptr++, remaining--)
  {
    uniquelength = substringinfo->uniqueforward(substringinfo->genericindex,
                                                &witnessposition,
                                                qptr,
                                                query+querylen);
    if (uniquelength > 0)
    {
      if (substringinfo->processuniquelength(substringinfo->alphabet,
                                             query,
                                             uniquelength,
                                             (unsigned long) (qptr-query),
                                             witnessposition,
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
                               Seqpos refpos,
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
    if (rangespecinfo->showrefpos)
    {
      printf(" " FormatSeqpos,PRINTSeqposcast(refpos));
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

int findsubqueryuniqueforward(const void *genericindex,
                              Uniqueforwardfunction uniqueforward,
                              const Alphabet *alphabet,
                              const StrArray *queryfilenames,
                              Definedunsignedlong minlength,
                              Definedunsignedlong maxlength,
                              bool showsequence,
                              bool showquerypos,
                              bool showrefpos,
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
  rangespecinfo.showrefpos = showrefpos;
  substringinfo.preprocessuniquelength = showunitnum;
  substringinfo.processuniquelength = showifinlengthrange;
  substringinfo.postprocessuniquelength = NULL;
  substringinfo.alphabet = alphabet;
  substringinfo.processinfo = &rangespecinfo;
  substringinfo.uniqueforward = uniqueforward;
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
  return haserr ? -1 : 0;
}
