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
#include <stdbool.h>
#include "core/alphabet.h"
#include "core/error.h"
#include "core/seq_iterator_sequence_buffer_api.h"
#include "core/unused_api.h"
#include "core/defined-types.h"
#include "core/codetype.h"
#include "core/encseq.h"
#include "core/format64.h"
#include "core/ma_api.h"
#include "optionargmode.h"
#include "greedyfwdmat.h"
#include "initbasepower.h"

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
typedef void (*Processgmatchlength)(const GtAlphabet *,
                                    const GtUchar *,
                                    unsigned long,
                                    unsigned long,
                                    unsigned long,
                                    void *);
typedef void (*Postprocessgmatchlength)(const GtAlphabet *,
                                        uint64_t,
                                        const char *,
                                        const GtUchar *,
                                        unsigned long,
                                        void *);

typedef struct
{
  const void *genericindex;
  unsigned long totallength;
  const GtAlphabet *alphabet;
  Greedygmatchforwardfunction gmatchforward;
  Preprocessgmatchlength preprocessgmatchlength;
  Processgmatchlength processgmatchlength;
  Postprocessgmatchlength postprocessgmatchlength;
  void *processinfo;
  const GtEncseq *encseq;
} Substringinfo;

#ifndef NDEBUG
static void checkifsequenceisthere(const GtEncseq *encseq,
                                   unsigned long witnessposition,
                                   unsigned long gmatchlength,
                                   const GtUchar *qptr)
{
  unsigned long i;
  GtUchar cc;

  for (i=0; i<gmatchlength; i++)
  {
    cc = gt_encseq_get_encoded_char_nospecial(encseq,
                                              witnessposition+i,
                                              GT_READMODE_FORWARD);
    if (qptr[i] != cc)
    {
      fprintf(stderr,"sequence of length %lu at witnesspos %lu"
                     " query[%lu] = %u != %u = subject[%lu]\n",
                     gmatchlength,
                     witnessposition,
                     i,
                     (unsigned int) qptr[i],
                     (unsigned int) cc,
                     witnessposition+(unsigned long) i);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
}
#endif

static void gmatchposinsinglesequence(Substringinfo *substringinfo,
                                      uint64_t unitnum,
                                      const GtUchar *query,
                                      unsigned long querylen,
                                      const char *desc)
{
  const GtUchar *qptr;
  unsigned long gmatchlength, remaining;
  unsigned long witnessposition, *wptr;

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
                                                0,
                                                0,
                                                substringinfo->totallength,
                                                wptr,
                                                qptr,
                                                query+querylen);
    if (gmatchlength > 0)
    {
#ifndef NDEBUG
      if (substringinfo->encseq != NULL)
      {
        gt_assert(wptr != NULL);
        checkifsequenceisthere(substringinfo->encseq,
                               witnessposition,
                               gmatchlength,
                               qptr);
      }
#endif
      substringinfo->processgmatchlength(substringinfo->alphabet,
                                         query,
                                         gmatchlength,
                                         (unsigned long) (qptr-query),
                                         wptr == NULL
                                           ? (unsigned long) 0
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
                        GT_UNUSED void *info)
{
  printf("unit " Formatuint64_t, PRINTuint64_tcast(unitnum));
  if (desc != NULL && desc[0] != '\0')
  {
    printf(" (%s)",desc);
  }
  printf("\n");
}

static void showifinlengthrange(const GtAlphabet *alphabet,
                                const GtUchar *start,
                                unsigned long gmatchlength,
                                unsigned long querystart,
                                unsigned long subjectpos,
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
      printf(" %lu",subjectpos);
    }
    if (rangespecinfo->showsequence)
    {
      (void) putchar(' ');
      gt_alphabet_decode_seq_to_fp(alphabet,stdout,start + querystart,
                                   gmatchlength);
    }
    (void) putchar('\n');
  }
}

int gt_findsubquerygmatchforward(const GtEncseq *encseq,
                              const void *genericindex,
                              unsigned long totallength,
                              Greedygmatchforwardfunction gmatchforward,
                              const GtAlphabet *alphabet,
                              const GtStrArray *queryfilenames,
                              Definedunsignedlong minlength,
                              Definedunsignedlong maxlength,
                              bool showsequence,
                              bool showquerypos,
                              bool showsubjectpos,
                              GtError *err)
{
  Substringinfo substringinfo;
  Rangespecinfo rangespecinfo;
  bool haserr = false;
  GtSeqIterator *seqit;
  const GtUchar *query;
  unsigned long querylen;
  char *desc = NULL;
  int retval;
  uint64_t unitnum;

  gt_error_check(err);
  substringinfo.genericindex = genericindex;
  substringinfo.totallength = totallength;
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
  seqit = gt_seq_iterator_sequence_buffer_new(queryfilenames, err);
  if (!seqit)
    haserr = true;
  if (!haserr)
  {
    gt_seq_iterator_set_symbolmap(seqit, gt_alphabet_symbolmap(alphabet));
    for (unitnum = 0; /* Nothing */; unitnum++)
    {
      retval = gt_seq_iterator_next(seqit,
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
    }
    gt_seq_iterator_delete(seqit);
  }
  return haserr ? -1 : 0;
}

#ifdef WITHrunsubstringiteration
int runsubstringiteration(Greedygmatchforwardfunction gmatchforward,
                          const void *genericindex,
                          unsigned long totalwidth,
                          const unsigned long *leftborder,
                          const unsigned long *countspecialcodes,
                          const Alphabet *alphabet,
                          unsigned int prefixlength,
                          const GtStrArray *queryfilenames,
                          GtError *err)
{
  Substriter *substriter;
  Substring substring;
  bool haserr = false;
  int retval;
  unsigned int numofchars;
  unsigned long gmatchlength, gmatchlength2;
  GtCodetype maxcode;
  GtBucketspecification bucketspec;

  substriter->seqit = gt_seq_iterator_new(filenames,
                                      getsymbolmapAlphabet(alphabet),
                                      true);
  substriter = gt_substriter_new(queryfilenames,alphabet,prefixlength);
  numofchars = getnumofcharsAlphabet(alphabet);
  maxcode = ontheflybasepower(numofchars,prefixlength);
  while (true)
  {
    retval = gt_substriter_next(&substring,substriter,err);
    if (retval < 0)
    {
      haserr = true;
      break;
    }
    if (retval == 0)
    {
      break;
    }
    gt_assert(substring.remaining >= (unsigned long) prefixlength);
    gmatchlength = gmatchforward(genericindex,
                                 0,
                                 0,
                                 totalwidth,
                                 NULL,
                                 substring.currentptr,
                                 substring.currentptr + substring.remaining);
    if (leftborder != NULL)
    {
      gt_bcktab_calcboundaries(&bucketspec,
                               leftborder,
                               countspecialcodes,
                               substring.currentcode,
                               maxcode,
                               totalwidth,
                               substring.currentcode % numofchars,
                               numofchars);
      if (bucketspec.nonspecialsinbucket > 0)
      {
        gmatchlength2 = gmatchforward(genericindex,
                                      (unsigned long) prefixlength,
                                      bucketspec.left,
                                      bucketspec.left
                                        + bucketspec.nonspecialsinbucket - 1,
                                      NULL,
                                      substring.currentptr+prefixlength,
                                      substring.currentptr+substring.remaining);
#ifndef NDEBUG
        if (gmatchlength2 != gmatchlength)
        {
          fprintf(stderr,"at offset %lu:\n",(unsigned long)
                                              (substring.currentptr -
                                               substring.start));
          fprintf(stderr,"bucketspec=(%lu,%lu)\n",
                          bucketspec.left,
                          bucketspec.left+bucketspec.nonspecialsinbucket-1);
          fprintf(stderr,"gmatchlength2 = %lu != %lu = gmatchlength\n",
                          gmatchlength2,gmatchlength);
          exit(GT_EXIT_PROGRAMMING_ERROR);
        }
#endif
      }
    }
  }
  gt_substriter_delete(substriter);
  return haserr ? -1 : 0;
}

int runsubstringiteration(Greedygmatchforwardfunction gmatchforward,
                          const void *genericindex,
                          unsigned long totalwidth,
                          const unsigned long *leftborder,
                          const unsigned long *countspecialcodes,
                          const Alphabet *alphabet,
                          unsigned int prefixlength,
                          const GtStrArray *queryfilenames,
                          GtError *err)
{
  GtSeqIterator *seqit;
  const GtUchar *query;
  unsigned long querylen;
  char *desc = NULL;
  Substriter *substriter;
  Substring substring;
  bool haserr = false;
  int retval;
  unsigned int numofchars;
  unsigned long gmatchlength, gmatchlength2;
  GtCodetype maxcode;
  GtBucketspecification bucketspec;
  bool haserr = false;

  seqit = gt_seq_iterator_new(queryfilenames,getsymbolmapAlphabet(alphabet),
  *                          true);
  for (unitnum = 0; ; unitnum++)
  {
    retval = gt_seq_iterator_next(seqit,
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
  }
  gt_seq_iterator_delete(seqit);
  return haserr ? -1 : 0;
}
#endif
