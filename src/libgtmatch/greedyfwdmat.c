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
#include "greedyfwdmat.h"
#include "substriter.h"
#include "encseq-def.h"

typedef struct
{
  bool showsequence,
       showquerypos,
       showsubjectpos;
  Definedunsignedlong minlength,
                      maxlength;
} Rangespecinfo;

typedef void (*Preprocessgmatchlength)(uint64_t,
                                       const char *,
                                       void *);
typedef void (*Processgmatchlength)(const Alphabet *,
                                    const Uchar *,
                                    unsigned long,
                                    unsigned long,
                                    Seqpos,
                                    void *);
typedef void (*Postprocessgmatchlength)(const Alphabet *,
                                        uint64_t,
                                        const char *,
                                        const Uchar *,
                                        unsigned long,
                                        void *);

typedef struct
{
  const void *genericindex;
  const Alphabet *alphabet;
  Greedygmatchforwardfunction gmatchforward;
  Preprocessgmatchlength preprocessgmatchlength;
  Processgmatchlength processgmatchlength;
  Postprocessgmatchlength postprocessgmatchlength;
  void *processinfo;
  const Encodedsequence *encseq;
} Substringinfo;

static void checkifsequenceisthere(const Encodedsequence *encseq,
                                   Seqpos witnessposition,
                                   unsigned long gmatchlength,
                                   const Uchar *qptr)
{
  unsigned long i;
  Uchar cc;

  for (i=0; i<gmatchlength; i++)
  {
    cc = getencodedchar(encseq,witnessposition+i,Forwardmode);
    if (qptr[i] != cc)
    {
      fprintf(stderr,"sequence of length %lu at witnesspos " FormatSeqpos
                     " query[%lu] = %u != %u = subject[" FormatSeqpos "]\n",
                     gmatchlength,
                     PRINTSeqposcast(witnessposition),
                     i,
                     (unsigned int) qptr[i],
                     (unsigned int) cc,
                     PRINTSeqposcast(witnessposition+(Seqpos) i));
      exit(EXIT_FAILURE); /* Program error */
    }
  }
}

static void gmatchposinsinglesequence(Substringinfo *substringinfo,
                                      uint64_t unitnum,
                                      const Uchar *query,
                                      unsigned long querylen,
                                      const char *desc)
{
  const Uchar *qptr;
  unsigned long gmatchlength, remaining;
  Seqpos witnessposition, *wptr;

  if (substringinfo->preprocessgmatchlength != NULL)
  {
    substringinfo->preprocessgmatchlength(unitnum,
                                          desc,
                                          substringinfo->processinfo);
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
    gmatchlength = substringinfo->gmatchforward(substringinfo->genericindex,
                                                wptr,
                                                qptr,
                                                query+querylen);
    if (gmatchlength > 0)
    {
      if (substringinfo->encseq != NULL)
      {
        assert(wptr != NULL);
        checkifsequenceisthere(substringinfo->encseq,
                               witnessposition,
                               gmatchlength,
                               qptr);
      }
      substringinfo->processgmatchlength(substringinfo->alphabet,
                                         query,
                                         gmatchlength,
                                         (unsigned long) (qptr-query),
                                         wptr == NULL
                                           ? (Seqpos) 0
                                           : witnessposition,
                                         substringinfo->processinfo);
    }
  }
  if (substringinfo->postprocessgmatchlength != NULL)
  {
    substringinfo->postprocessgmatchlength(substringinfo->alphabet,
                                           unitnum,
                                           desc,
                                           query,
                                           querylen,
                                           substringinfo->processinfo);
  }
}

static void showunitnum(uint64_t unitnum,
                       const char *desc,
                       /*@unused@*/ void *info)
{
  printf("unit " Formatuint64_t, PRINTuint64_tcast(unitnum));
  if (desc != NULL && desc[0] != '\0')
  {
    printf(" (%s)",desc);
  }
  printf("\n");
}

static void showifinlengthrange(const Alphabet *alphabet,
                                const Uchar *start,
                                unsigned long gmatchlength,
                                unsigned long querystart,
                                Seqpos subjectpos,
                                void *info)
{
  Rangespecinfo *rangespecinfo = (Rangespecinfo *) info;

  if ((!rangespecinfo->minlength.defined ||
      gmatchlength >= rangespecinfo->minlength.valueunsignedlong) &&
     (!rangespecinfo->maxlength.defined ||
      gmatchlength <= rangespecinfo->maxlength.valueunsignedlong))
  {
    if (rangespecinfo->showquerypos)
    {
      printf("%lu ",querystart);
    }
    printf("%lu",gmatchlength);
    if (rangespecinfo->showsubjectpos)
    {
      printf(" " FormatSeqpos,PRINTSeqposcast(subjectpos));
    }
    if (rangespecinfo->showsequence)
    {
      (void) putchar(' ');
      showsymbolstring(alphabet,
                       start + querystart,
                       gmatchlength);
    }
    (void) putchar('\n');
  }
}

int findsubquerygmatchforward(const Encodedsequence *encseq,
                              const void *genericindex,
                              Greedygmatchforwardfunction gmatchforward,
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
  substringinfo.preprocessgmatchlength = showunitnum;
  substringinfo.processgmatchlength = showifinlengthrange;
  substringinfo.postprocessgmatchlength = NULL;
  substringinfo.alphabet = alphabet;
  substringinfo.processinfo = &rangespecinfo;
  substringinfo.gmatchforward = gmatchforward;
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
    gmatchposinsinglesequence(&substringinfo,
                              unitnum,
                              query,
                              querylen,
                              desc);
    FREESPACE(desc);
  }
  seqiterator_delete(seqit);
  return haserr ? -1 : 0;
}

int runsubstringiteration(const Alphabet *alphabet,
                          const StrArray *queryfilenames,
                          unsigned int prefixlength,
                          Error *err)
{
  Substriter *substriter;
  Substring substring;
  bool haserr = false;
  int retval;

  substriter = substriter_new(queryfilenames,alphabet,prefixlength);
  while (true)
  {
    retval = substriter_next(&substring,substriter,err);
    if (retval < 0)
    {
      haserr = true;
      break;
    }
    if (retval == 0)
    {
      break;
    }
  }
  substriter_delete(&substriter);
  return haserr ? -1 : 0;
}
