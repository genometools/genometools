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

#include "core/chardef.h"
#include "core/divmodmul.h"
#include "core/minmax.h"
#include "core/seq_iterator_sequence_buffer_api.h"
#include "core/types_api.h"
#include "core/timer_api.h"
#include "core/format64.h"
#include "sarr-def.h"
#include "revcompl.h"
#include "lcpinterval.h"
#include "esa-map.h"
#include "echoseq.h"
#include "sfx-apfxlen.h"
#include "sfx-suffixer.h"
#include "esa-minunique.h"
#include "esa-mmsearch.h"

typedef struct
{
  const GtUchar *sequence;
  bool reversecopy;
  const GtEncseq *encseq;
  GtReadmode readmode;
  unsigned long startpos, length;
} GtQueryrep;

typedef struct
{
  const GtQueryrep *queryrep;
  unsigned long offset; /* position relative to startpos */
} GtQuerysubstring;

static GtUchar gt_mmsearch_accessquery(const GtQueryrep *queryrep,
                                       unsigned long pos)
{
  unsigned long abspos;

  gt_assert(queryrep != NULL && pos < queryrep->length);
  abspos = queryrep->startpos + pos;
  if (queryrep->sequence != NULL)
  {
    gt_assert(queryrep->readmode == GT_READMODE_FORWARD);
    return queryrep->sequence[abspos];
  } else
  {
    gt_assert(queryrep->readmode != GT_READMODE_FORWARD &&
              queryrep->encseq != NULL);
    return gt_encseq_get_encoded_char(queryrep->encseq,abspos,
                                      queryrep->readmode);
  }
}

#define GT_MMSEARCH_COMPARE(OFFSET,LCPLEN)\
        sidx = (OFFSET) + (LCPLEN);\
        if (sidx < totallength)\
        {\
          gt_encseq_reader_reinit_with_readmode(esr,dbencseq,readmode,sidx);\
        }\
        for (/* Nothing */ ; /* Nothing */; sidx++, (LCPLEN)++)\
        {\
          if ((LCPLEN) >= (unsigned long) minmatchlength)\
          {\
            retcode = 0;\
            break;\
          }\
          if (sidx >= totallength)\
          {\
            retcode = -1;\
            break;\
          }\
          currentdbchar = gt_encseq_reader_next_encoded_char(esr);\
          currentquerychar = gt_mmsearch_accessquery(querysubstring->queryrep,\
                                         querysubstring->offset + (LCPLEN));\
          retcode = (int) (currentquerychar - currentdbchar);\
          if (retcode == 0)\
          {\
            if (ISSPECIAL(currentdbchar) && ISSPECIAL(currentquerychar))\
            {\
              retcode = -1;\
              break;\
            }\
          } else\
          {\
            break;\
          }\
        }

static bool gt_mmsearch(const GtEncseq *dbencseq,
                        GtEncseqReader *esr,
                        const ESASuffixptr *suftab,
                        GtReadmode readmode,
                        Lcpinterval *lcpitv,
                        const GtQuerysubstring *querysubstring,
                        unsigned long minmatchlength)
{
  unsigned long left, leftsave, mid, right, lpref, rpref,
                totallength, lcplen, sidx;
  int retcode = 0;
  GtUchar currentdbchar, currentquerychar;

  totallength = gt_encseq_total_length(dbencseq);
  leftsave = left = lcpitv->left;
  right = lcpitv->right;
  lcplen = lcpitv->offset;
  GT_MMSEARCH_COMPARE(ESASUFFIXPTRGET(suftab,left),lcplen);
  if (retcode > 0)
  {
    lpref = lcplen;
    lcplen = lcpitv->offset;
    GT_MMSEARCH_COMPARE(ESASUFFIXPTRGET(suftab,right),lcplen);
    if (retcode > 0)
    {
      return false;
    } else
    {
      rpref = lcplen;
      while (right > left + 1)
      {
        mid = GT_DIV2(left+right);
        lcplen = MIN(lpref,rpref);
        GT_MMSEARCH_COMPARE(ESASUFFIXPTRGET(suftab,mid),lcplen);
        if (retcode <= 0)
        {
          right = mid;
          rpref = lcplen;
        } else
        {
          left = mid;
          lpref = lcplen;
        }
      }
      lcpitv->left = right;
    }
  }

  left = leftsave;
  right = lcpitv->right;
  lcplen = lcpitv->offset;
  GT_MMSEARCH_COMPARE(ESASUFFIXPTRGET(suftab,left),lcplen);
  if (retcode < 0)
  {
    return false;
  } else
  {
    lpref = lcplen;
    lcplen = lcpitv->offset;
    GT_MMSEARCH_COMPARE(ESASUFFIXPTRGET(suftab,right),lcplen);
    if (retcode >= 0)
    {
      lcpitv->right = right;
    } else
    {
      rpref = lcplen;
      while (right > left + 1)
      {
        mid = GT_DIV2(left+right);
        lcplen = MIN(lpref,rpref);
        GT_MMSEARCH_COMPARE(ESASUFFIXPTRGET(suftab,mid),lcplen);
        if (retcode >= 0)
        {
          left = mid;
          lpref = lcplen;
        } else
        {
          right = mid;
          rpref = lcplen;
        }
      }
      lcpitv->right = left;
    }
  }
  return true;
}

struct GtMMsearchiterator
{
  Lcpinterval lcpitv;
  unsigned long sufindex;
  const ESASuffixptr *suftab;
  GtEncseqReader *esr;
};

static GtMMsearchiterator *gt_mmsearchiterator_new_generic(
                                       const GtEncseq *dbencseq,
                                       const ESASuffixptr *suftab,
                                       unsigned long leftbound,
                                       unsigned long rightbound,
                                       unsigned long itvoffset,
                                       GtReadmode readmode,
                                       const GtQuerysubstring *querysubstring,
                                       unsigned long minmatchlength)

{
  GtMMsearchiterator *mmsi = gt_malloc(sizeof *mmsi);

  mmsi->lcpitv.left = leftbound;
  mmsi->lcpitv.right = rightbound;
  mmsi->lcpitv.offset = itvoffset;
  mmsi->suftab = suftab;
  mmsi->esr = gt_encseq_create_reader_with_readmode(dbencseq, readmode, 0);
  if (!gt_mmsearch(dbencseq,mmsi->esr,suftab,readmode,&mmsi->lcpitv,
                   querysubstring,minmatchlength))
  {
    mmsi->lcpitv.left = 1UL;
    mmsi->lcpitv.right = 0;
  }
  mmsi->sufindex = mmsi->lcpitv.left;
  return mmsi;
}

GtMMsearchiterator *gt_mmsearchiterator_new_complete_olain(
                                   const GtEncseq *dbencseq,
                                   const void *voidsuftab, /* XXX */
                                   unsigned long leftbound,
                                   unsigned long rightbound,
                                   unsigned long itvoffset,
                                   GtReadmode readmode,
                                   const GtUchar *pattern,
                                   unsigned long patternlength)
{
  GtQueryrep queryrep;
  GtQuerysubstring querysubstring;
  const ESASuffixptr *suftab = (const ESASuffixptr *) voidsuftab; /* XXX */

  queryrep.sequence = pattern;
  queryrep.reversecopy = false;
  queryrep.encseq = NULL;
  queryrep.readmode = GT_READMODE_FORWARD;
  queryrep.startpos = 0;
  queryrep.length = (unsigned long) patternlength;
  querysubstring.queryrep = &queryrep;
  querysubstring.offset = 0;
  return gt_mmsearchiterator_new_generic(dbencseq,
                                         suftab,
                                         leftbound,
                                         rightbound,
                                         itvoffset,
                                         readmode,
                                         &querysubstring,
                                         patternlength);
}

unsigned long gt_mmsearchiterator_count(const GtMMsearchiterator *mmsi)
{
  if (mmsi->lcpitv.left > mmsi->lcpitv.right)
  {
    return 0;
  }
  return mmsi->lcpitv.right - mmsi->lcpitv.left + 1;
}

bool gt_mmsearchiterator_next(unsigned long *dbstart,GtMMsearchiterator *mmsi)
{
  if (mmsi->sufindex <= mmsi->lcpitv.right)
  {
    *dbstart = ESASUFFIXPTRGET(mmsi->suftab,mmsi->sufindex++);
    return true;
  }
  return false;
}

bool gt_mmsearchiterator_isempty(const GtMMsearchiterator *mmsi)
{
  return mmsi == NULL || mmsi->lcpitv.left > mmsi->lcpitv.right;
}

bool gt_mmsearchiterator_identical(const GtMMsearchiterator *mmsi1,
                                   const GtMMsearchiterator *mmsi2)
{
  gt_assert(mmsi1 != NULL);
  gt_assert(mmsi2 != NULL);
  return mmsi1->lcpitv.left == mmsi2->lcpitv.left &&
         mmsi1->lcpitv.right == mmsi2->lcpitv.right;
}

void gt_mmsearchiterator_delete(GtMMsearchiterator *mmsi)
{
  if (mmsi != NULL)
  {
    gt_encseq_reader_delete(mmsi->esr);
    gt_free(mmsi);
  }
}

static bool gt_mmsearch_isleftmaximal(const GtEncseq *dbencseq,
                                      GtReadmode readmode,
                                      unsigned long dbstart,
                                      const GtQuerysubstring *querysubstring)
{
  GtUchar dbleftchar;

  if (dbstart == 0 || querysubstring->offset == 0)
  {
    return true;
  }
  dbleftchar = gt_encseq_get_encoded_char(dbencseq, /* Random access */
                              dbstart-1,
                              readmode);
  if (ISSPECIAL(dbleftchar) ||
      dbleftchar != gt_mmsearch_accessquery(querysubstring->queryrep,
                                            querysubstring->offset-1))
  {
    return true;
  }
  return false;
}

static bool gt_mum_isleftmaximal(const GtEncseq *dbencseq,
                                 GtReadmode readmode,
                                 unsigned long dbstart,
                                 unsigned long queryoffset,
                                 const GtUchar *query)
{
  GtUchar dbleftchar;

  if (dbstart == 0 || queryoffset == 0)
  {
    return true;
  }
  dbleftchar = gt_encseq_get_encoded_char(dbencseq, /* Random access */
                                          dbstart-1,
                                          readmode);
  if (ISSPECIAL(dbleftchar) || dbleftchar != query[queryoffset-1])
  {
    return true;
  }
  return false;
}

static unsigned long gt_mmsearch_extendright(const GtEncseq *dbencseq,
                                             GtEncseqReader *esr,
                                             GtReadmode readmode,
                                             unsigned long totallength,
                                             unsigned long dbend,
                                             const GtQuerysubstring
                                               *querysubstring,
                                             unsigned long matchlength)
{
  GtUchar dbchar;
  unsigned long dbpos, querypos;

  if (dbend < totallength)
  {
    gt_encseq_reader_reinit_with_readmode(esr,dbencseq,readmode,dbend);
  }
  for (dbpos = dbend, querypos = querysubstring->offset + matchlength;
       dbpos < totallength &&
       querypos < querysubstring->queryrep->length;
       dbpos++, querypos++)
  {
    dbchar = gt_encseq_reader_next_encoded_char(esr);
    if (ISSPECIAL(dbchar) ||
        dbchar != gt_mmsearch_accessquery(querysubstring->queryrep,querypos))
    {
      break;
    }
  }
  return dbpos - dbend;
}

static int gt_queryuniquematch(bool selfmatch,
                               const Suffixarray *suffixarray,
                               uint64_t queryunitnum,
                               const GtQueryrep *queryrep,
                               unsigned long minmatchlength,
                               GtProcessquerymatch processquerymatch,
                               void *processquerymatchinfo,
                               GtQuerymatch *querymatchspaceptr,
                               GtError *err)
{
  unsigned long offset,
                totallength = gt_encseq_total_length(suffixarray->encseq),
                localqueryoffset = 0;
  uint64_t localqueryunitnum = queryunitnum;
  bool haserr = false;

  gt_assert(!selfmatch && queryrep->length >= minmatchlength);
  for (offset = 0; offset <= queryrep->length - minmatchlength; offset++)
  {
    unsigned long matchlen, dbstart;

    matchlen = gt_suffixarrayfindmums (suffixarray,
                                       0,
                                       0, /* leftbound */
                                       totallength, /* rightbound */
                                       &dbstart,
                                       queryrep->sequence + offset,
                                       queryrep->sequence + queryrep->length);
    if (dbstart != ULONG_MAX &&
        matchlen >= minmatchlength &&
        gt_mum_isleftmaximal(suffixarray->encseq,
                             suffixarray->readmode,
                             dbstart,
                             offset,
                             queryrep->sequence))
    {
      gt_querymatch_fill(querymatchspaceptr,
                         matchlen,
                         dbstart,
                         queryrep->readmode,
                         queryrep->reversecopy,
                         0, /* score */
                         0, /* edist */
                         selfmatch,
                         localqueryunitnum,
                         matchlen,
                         localqueryoffset);
      if (processquerymatch(processquerymatchinfo,
                            suffixarray->encseq,
                            querymatchspaceptr,
                            queryrep->sequence,
                            queryrep->length,
                            err) != 0)
      {
        haserr = true;
      }
    }
    if (!haserr)
    {
      if (queryrep->sequence[offset] == (GtUchar) SEPARATOR)
      {
        localqueryunitnum++;
        localqueryoffset = 0;
      } else
      {
        localqueryoffset++;
      }
    }
  }
  return haserr ? -1 : 0;
}

static int gt_querysubstringmatch_generic(
                                     bool selfmatch,
                                     const GtEncseq *dbencseq,
                                     const ESASuffixptr *suftabpart,
                                     GtReadmode readmode,
                                     unsigned long numberofsuffixes,
                                     uint64_t queryunitnum,
                                     const GtQueryrep *queryrep,
                                     unsigned long minmatchlength,
                                     GtProcessquerymatch processquerymatch,
                                     void *processquerymatchinfo,
                                     GtQuerymatch *querymatchspaceptr,
                                     GtError *err)
{
  GtMMsearchiterator *mmsi;
  unsigned long totallength, localqueryoffset = 0;
  uint64_t localqueryunitnum = queryunitnum;
  GtQuerysubstring querysubstring;
  bool haserr = false;

  gt_assert(numberofsuffixes > 0);
  totallength = gt_encseq_total_length(dbencseq);
  querysubstring.queryrep = queryrep;
  for (querysubstring.offset = 0;
       querysubstring.offset <= queryrep->length - minmatchlength;
       querysubstring.offset++)
  {
    unsigned long dbstart;

    mmsi = gt_mmsearchiterator_new_generic(dbencseq,
                                           suftabpart,
                                           0, /* leftbound */
                                           numberofsuffixes-1, /* rightbound */
                                           0, /* offset */
                                           readmode,
                                           &querysubstring,
                                           minmatchlength);
    while (!haserr && gt_mmsearchiterator_next(&dbstart,mmsi))
    {
      if (gt_mmsearch_isleftmaximal(dbencseq,
                                    readmode,
                                    dbstart,
                                    &querysubstring))
      {
        unsigned long extend = gt_mmsearch_extendright(dbencseq,
                                                       mmsi->esr,
                                                       readmode,
                                                       totallength,
                                                       dbstart + minmatchlength,
                                                       &querysubstring,
                                                       minmatchlength);
        gt_querymatch_fill(querymatchspaceptr,
                           minmatchlength + extend,
                           dbstart,
                           queryrep->readmode,
                           queryrep->reversecopy,
                           0, /* score */
                           0, /* edist */
                           selfmatch,
                           localqueryunitnum,
                           minmatchlength + extend,
                           localqueryoffset);
        if (processquerymatch(processquerymatchinfo,
                              dbencseq,
                              querymatchspaceptr,
                              queryrep->sequence,
                              queryrep->length,
                              err) != 0)
        {
          haserr = true;
        }
      }
    }
    gt_mmsearchiterator_delete(mmsi);
    mmsi = NULL;
    if (!haserr)
    {
      if (gt_mmsearch_accessquery(queryrep,querysubstring.offset)
          == (GtUchar) SEPARATOR)
      {
        localqueryunitnum++;
        localqueryoffset = 0;
      } else
      {
        localqueryoffset++;
      }
    }
  }
  return haserr ? -1 : 0;
}

static int gt_querysubstringmatch(bool selfmatch,
                                  const Suffixarray *suffixarray,
                                  uint64_t queryunitnum,
                                  const GtQueryrep *queryrep,
                                  unsigned long minmatchlength,
                                  GtProcessquerymatch processquerymatch,
                                  void *processquerymatchinfo,
                                  GtQuerymatch *querymatchspaceptr,
                                  GtError *err)
{
  return gt_querysubstringmatch_generic(selfmatch,
                              suffixarray->encseq,
                              suffixarray->suftab,
                              suffixarray->readmode,
                              gt_encseq_total_length(suffixarray->encseq) + 1,
                              queryunitnum,
                              queryrep,
                              minmatchlength,
                              processquerymatch,
                              processquerymatchinfo,
                              querymatchspaceptr,
                              err);
}

typedef int (*GtQuerysubstringmatchfunc)(bool,
                                         const Suffixarray *,
                                         uint64_t,
                                         const GtQueryrep *,
                                         unsigned long,
                                         GtProcessquerymatch,
                                         void *,
                                         GtQuerymatch *,
                                         GtError *);

static int gt_callenumquerymatches_withindex(
                            GtQuerysubstringmatchfunc findquerymatches,
                            const Suffixarray *suffixarray,
                            const GtStrArray *queryfiles,
                            bool forwardstrand,
                            bool reversestrand,
                            unsigned int userdefinedleastlength,
                            GtProcessquerybeforematching
                               processquerybeforematching,
                            GtProcessquerymatch processquerymatch,
                            void *processquerymatchinfo,
                            GtError *err)
{
  GtSeqIterator *seqit;
  bool haserr = false;

  seqit = gt_seq_iterator_sequence_buffer_new(queryfiles, err);
  if (seqit == NULL)
  {
    haserr = true;
  } else
  {
    GtQuerymatch *querymatchspaceptr = gt_querymatch_new();
    const GtUchar *query;
    unsigned long querylen;
    int retval;
    uint64_t queryunitnum;
    GtUchar *queryreverse = NULL;
    unsigned long queryreverse_length = 0;
    char *desc = NULL;
    int mode;

    gt_seq_iterator_set_symbolmap(seqit,
                    gt_alphabet_symbolmap(gt_encseq_alphabet(
                                                        suffixarray->encseq)));
    for (queryunitnum = 0; /* Nothing */; queryunitnum++)
    {
      retval = gt_seq_iterator_next(seqit, &query, &querylen, &desc, err);
      if (retval < 0)
      {
        haserr = true;
        break;
      }
      if (retval == 0)
      {
        break;
      }
      if (querylen >= (unsigned long) userdefinedleastlength)
      {
        GtQueryrep queryrep;

        queryrep.encseq = NULL;
        queryrep.readmode = GT_READMODE_FORWARD;
        queryrep.startpos = 0;
        queryrep.length = querylen;
        for (mode = 0; mode <= 1; mode++)
        {
          if (mode == 0 && forwardstrand)
          {
            queryrep.sequence = query;
            queryrep.reversecopy = false;
            if (processquerybeforematching != NULL)
            {
              processquerybeforematching(processquerymatchinfo,desc,query,
                                         querylen,true);
            }
          } else
          {
            if (mode == 1 && reversestrand)
            {
              if (querylen > queryreverse_length)
              {
                queryreverse = gt_realloc(queryreverse,
                                          sizeof (*queryreverse) * querylen);
                queryreverse_length = querylen;
              }
              gt_copy_reversecomplement(queryreverse,query,querylen);
              queryrep.sequence = queryreverse;
              queryrep.reversecopy = true;
              if (processquerybeforematching != NULL)
              {
                processquerybeforematching(processquerymatchinfo,desc,
                                           queryreverse,querylen,false);
              }
            } else
            {
              queryrep.sequence = NULL;
              queryrep.reversecopy = false;
            }
          }
          if (queryrep.sequence != NULL)
          {
            int ret = findquerymatches(false,
                                       suffixarray,
                                       queryunitnum,
                                       &queryrep,
                                       (unsigned long) userdefinedleastlength,
                                       processquerymatch,
                                       processquerymatchinfo,
                                       querymatchspaceptr,
                                       err);
            if (ret != 0)
            {
              haserr = true;
              break;
            }
          }
        }
      }
    }
    gt_seq_iterator_delete(seqit);
    gt_free(queryreverse);
    gt_querymatch_delete(querymatchspaceptr);
  }
  return haserr ? -1 : 0;
}

int gt_callenumquerymatches(const char *indexname,
                            const GtStrArray *queryfiles,
                            bool findmums,
                            bool forwardstrand,
                            bool reversestrand,
                            unsigned int userdefinedleastlength,
                            GtProcessquerybeforematching
                              processquerybeforematching,
                            GtProcessquerymatch processquerymatch,
                            void *processquerymatchinfo,
                            GtLogger *logger,
                            GtError *err)
{
  Suffixarray suffixarray;
  bool haserr = false;

  if (gt_mapsuffixarray(&suffixarray,
                        SARR_ESQTAB | SARR_SUFTAB | SARR_SSPTAB,
                        indexname,
                        logger,
                        err) != 0)
  {
    haserr = true;
  } else
  {
    if (gt_callenumquerymatches_withindex(findmums
                                            ? gt_queryuniquematch
                                            : gt_querysubstringmatch,
                                          &suffixarray,
                                          queryfiles,
                                          forwardstrand,
                                          reversestrand,
                                          userdefinedleastlength,
                                          processquerybeforematching,
                                          processquerymatch,
                                          processquerymatchinfo,
                                          err) != 0)
    {
      haserr = true;
    }
  }
  gt_freesuffixarray(&suffixarray);
  return haserr ? -1 : 0;
}

int gt_callenumselfmatches(const char *indexname,
                           GtReadmode queryreadmode,
                           unsigned int userdefinedleastlength,
                           GtProcessquerymatch processquerymatch,
                           void *processquerymatchinfo,
                           GtLogger *logger,
                           GtError *err)
{
  Suffixarray suffixarray;
  bool haserr = false;

  gt_assert(queryreadmode != GT_READMODE_FORWARD);
  if (gt_mapsuffixarray(&suffixarray,
                        SARR_ESQTAB | SARR_SUFTAB | SARR_SSPTAB,
                        indexname,
                        logger,
                        err) != 0)
  {
    haserr = true;
  } else
  {
    unsigned long seqnum, numofsequences, seqlength, seqstartpos;
    GtQuerymatch *querymatchspaceptr = gt_querymatch_new();
    GtQueryrep queryrep;

    numofsequences = gt_encseq_num_of_sequences(suffixarray.encseq);
    queryrep.sequence = NULL;
    queryrep.reversecopy = false;
    queryrep.encseq = suffixarray.encseq;
    queryrep.readmode = queryreadmode;
    for (seqnum = 0; seqnum < numofsequences; seqnum++)
    {
      seqstartpos = gt_encseq_seqstartpos(suffixarray.encseq, seqnum);
      seqlength = gt_encseq_seqlength(suffixarray.encseq, seqnum);
      if (seqlength >= (unsigned long) userdefinedleastlength)
      {
        queryrep.startpos = seqstartpos;
        queryrep.length = seqlength;
        if (gt_querysubstringmatch(true,
                                   &suffixarray,
                                   (uint64_t) seqnum,
                                   &queryrep,
                                   (unsigned long) userdefinedleastlength,
                                   processquerymatch,
                                   processquerymatchinfo,
                                   querymatchspaceptr,
                                   err) != 0)
        {
          haserr = true;
          break;
        }
      }
    }
    gt_querymatch_delete(querymatchspaceptr);
  }
  gt_freesuffixarray(&suffixarray);
  return haserr ? -1 : 0;
}

static int gt_constructsarrandrunmmsearch(
                 const GtEncseq *dbencseq,
                 GtReadmode readmode,
                 unsigned int prefixlength,
                 unsigned int numofparts,
                 unsigned long maximumspace,
                 const GtUchar *query,
                 unsigned long querylen,
                 unsigned int minlength,
                 GtProcessquerymatch processquerymatch,
                 void *processquerymatchinfo,
                 GtTimer *sfxprogress,
                 bool withprogressbar,
                 GtError *err)
{
  bool haserr = false;
  Sfxiterator *sfi;
  Sfxstrategy sfxstrategy;

  defaultsfxstrategy(&sfxstrategy,
                     gt_encseq_bitwise_cmp_ok(dbencseq) ? false : true);
  sfxstrategy.outsuftabonfile = false;
  sfi = gt_Sfxiterator_new(dbencseq,
                           readmode,
                           prefixlength,
                           numofparts,
                           maximumspace,
                           &sfxstrategy, /* sfxstrategy */
                           sfxprogress,
                           withprogressbar,
                           NULL, /* logger */
                           err);
  if (sfi == NULL)
  {
    haserr = true;
  } else
  {
    const GtSuffixsortspace *suffixsortspace;
    unsigned long numberofsuffixes;
    GtQuerymatch *querymatchspaceptr = gt_querymatch_new();
    GtQueryrep queryrep;

    queryrep.sequence = query;
    queryrep.reversecopy = false;
    queryrep.encseq = NULL;
    queryrep.readmode = GT_READMODE_FORWARD;
    queryrep.startpos = 0;
    queryrep.length = (unsigned long) querylen;
    while (true)
    {
      suffixsortspace = gt_Sfxiterator_next(&numberofsuffixes,NULL,sfi);
      if (suffixsortspace == NULL)
      {
        break;
      }
      if (gt_querysubstringmatch_generic(
                                false,
                                dbencseq,
                                (const ESASuffixptr *)
                                gt_suffixsortspace_ulong_get(suffixsortspace),
                                readmode,
                                numberofsuffixes,
                                0,
                                &queryrep,
                                (unsigned long) minlength,
                                processquerymatch,
                                processquerymatchinfo,
                                querymatchspaceptr,
                                err) != 0)
      {
        haserr = true;
        break;
      }
    }
    gt_querymatch_delete(querymatchspaceptr);
  }
  if (gt_Sfxiterator_delete(sfi,err) != 0)
  {
    haserr = true;
  }
  return haserr ? -1 : 0;
}

int gt_sarrquerysubstringmatch(const GtUchar *dbseq,
                               unsigned long dblen,
                               const GtUchar *query,
                               unsigned long querylen,
                               unsigned int minlength,
                               GtAlphabet *alpha,
                               GtProcessquerymatch processquerymatch,
                               void *processquerymatchinfo,
                               GtLogger *logger,
                               GtError *err)
{
  unsigned int numofchars, recommendedprefixlength;
  bool haserr = false;
  GtEncseq *dbencseq;
  GtEncseqBuilder *eb;

  eb = gt_encseq_builder_new(alpha);
  gt_encseq_builder_disable_multiseq_support(eb);
  gt_encseq_builder_disable_description_support(eb);
  gt_encseq_builder_set_logger(eb, logger);
  gt_encseq_builder_add_encoded(eb, dbseq, dblen, NULL);
  dbencseq = gt_encseq_builder_build(eb, err);
  gt_encseq_builder_delete(eb);

  numofchars = gt_alphabet_num_of_chars(alpha);
  recommendedprefixlength
    = gt_recommendedprefixlength(numofchars,dblen,
                                 GT_RECOMMENDED_MULTIPLIER_DEFAULT,
                                 true);
  if (gt_constructsarrandrunmmsearch(dbencseq,
                                     GT_READMODE_FORWARD,
                                     recommendedprefixlength,
                                     1U, /* parts */
                                     0, /* maximumspace */
                                     query,
                                     querylen,
                                     minlength,
                                     processquerymatch,
                                     processquerymatchinfo,
                                     NULL,
                                     false,
                                     err) != 0)
  {
    haserr = true;
  }
  gt_encseq_delete(dbencseq);
  dbencseq = NULL;
  return haserr ? -1 : 0;
}
