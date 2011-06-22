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
#include "core/seqiterator_sequence_buffer.h"
#include "core/types_api.h"

#include "sarr-def.h"
#include "lcpinterval.h"
#include "spacedef.h"
#include "esa-mmsearch.h"
#include "core/timer_api.h"
#include "core/format64.h"
#include "stamp.h"
#include "esa-map.h"
#include "echoseq.h"
#include "sfx-apfxlen.h"
#include "sfx-suffixer.h"

typedef struct
{
  const GtUchar *sequence;
  const GtEncseq *encseq;
  GtReadmode readmode;
  unsigned long startpos, length;
} Queryrep;

typedef struct
{
  const Queryrep *queryrep;
  unsigned long offset; /* position relative to startpos */
} Querysubstring;

static GtUchar accessquery(const Queryrep *queryrep,unsigned long pos)
{
  unsigned long abspos = queryrep->startpos + pos;

  gt_assert(pos < queryrep->length);
  if (queryrep->sequence != NULL)
  {
    gt_assert(queryrep->readmode == GT_READMODE_FORWARD);
    return queryrep->sequence[abspos];
  } else
  {
    gt_assert(queryrep->readmode != GT_READMODE_FORWARD);
    gt_assert(queryrep->encseq != NULL);
    return gt_encseq_get_encoded_char(queryrep->encseq,
                                             abspos,
                                             queryrep->readmode);
  }
}

#define COMPARE(OFFSET,LCPLEN)\
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
          currentquerychar = accessquery(querysubstring->queryrep,\
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

static bool mmsearch(const GtEncseq *dbencseq,
                     GtEncseqReader *esr,
                     const ESASuffixptr *suftab,
                     GtReadmode readmode,
                     Lcpinterval *lcpitv,
                     const Querysubstring *querysubstring,
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
  COMPARE(ESASUFFIXPTRGET(suftab,left),lcplen);
  if (retcode > 0)
  {
    lpref = lcplen;
    lcplen = lcpitv->offset;
    COMPARE(ESASUFFIXPTRGET(suftab,right),lcplen);
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
        COMPARE(ESASUFFIXPTRGET(suftab,mid),lcplen);
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
  COMPARE(ESASUFFIXPTRGET(suftab,left),lcplen);
  if (retcode < 0)
  {
    return false;
  } else
  {
    lpref = lcplen;
    lcplen = lcpitv->offset;
    COMPARE(ESASUFFIXPTRGET(suftab,right),lcplen);
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
        COMPARE(ESASUFFIXPTRGET(suftab,mid),lcplen);
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

 struct MMsearchiterator
{
  Lcpinterval lcpitv;
  unsigned long sufindex;
  const ESASuffixptr *suftab;
  GtEncseqReader *esr;
};

static MMsearchiterator *newmmsearchiterator_generic(
                                       const GtEncseq *dbencseq,
                                       const ESASuffixptr *suftab,
                                       unsigned long leftbound,
                                       unsigned long rightbound,
                                       unsigned long itvoffset,
                                       GtReadmode readmode,
                                       const Querysubstring *querysubstring,
                                       unsigned long minmatchlength)

{
  MMsearchiterator *mmsi;

  ALLOCASSIGNSPACE(mmsi,NULL,MMsearchiterator,1);
  mmsi->lcpitv.left = leftbound;
  mmsi->lcpitv.right = rightbound;
  mmsi->lcpitv.offset = itvoffset;
  mmsi->suftab = suftab;
  mmsi->esr = gt_encseq_create_reader_with_readmode(dbencseq, readmode, 0);
  if (!mmsearch(dbencseq,mmsi->esr,suftab,readmode,&mmsi->lcpitv,
                querysubstring,minmatchlength))
  {
    mmsi->lcpitv.left = (unsigned long) 1;
    mmsi->lcpitv.right = 0;
  }
  mmsi->sufindex = mmsi->lcpitv.left;
  return mmsi;
}

MMsearchiterator *gt_newmmsearchiteratorcomplete_plain(
                                   const GtEncseq *dbencseq,
                                   const void *voidsuftab, /* XXX */
                                   unsigned long leftbound,
                                   unsigned long rightbound,
                                   unsigned long itvoffset,
                                   GtReadmode readmode,
                                   const GtUchar *pattern,
                                   unsigned long patternlength)
{
  Queryrep queryrep;
  Querysubstring querysubstring;
  const ESASuffixptr *suftab = (const ESASuffixptr *) voidsuftab; /* XXX */

  queryrep.sequence = pattern;
  queryrep.encseq = NULL;
  queryrep.readmode = GT_READMODE_FORWARD;
  queryrep.startpos = 0;
  queryrep.length = (unsigned long) patternlength;
  querysubstring.queryrep = &queryrep;
  querysubstring.offset = 0;
  return newmmsearchiterator_generic(dbencseq,
                                     suftab,
                                     leftbound,
                                     rightbound,
                                     itvoffset,
                                     readmode,
                                     &querysubstring,
                                     patternlength);
}

unsigned long gt_countmmsearchiterator(const MMsearchiterator *mmsi)
{
  if (mmsi->lcpitv.left > mmsi->lcpitv.right)
  {
    return 0;
  }
  return mmsi->lcpitv.right - mmsi->lcpitv.left + 1;
}

bool gt_nextmmsearchiterator(unsigned long *dbstart,MMsearchiterator *mmsi)
{
  if (mmsi->sufindex <= mmsi->lcpitv.right)
  {
    *dbstart = ESASUFFIXPTRGET(mmsi->suftab,mmsi->sufindex++);
    return true;
  }
  return false;
}

bool gt_isemptymmsearchiterator(const MMsearchiterator *mmsi)
{
  return mmsi == NULL || mmsi->lcpitv.left > mmsi->lcpitv.right;
}

bool gt_identicalmmsearchiterators(const MMsearchiterator *mmsi1,
                                const MMsearchiterator *mmsi2)
{
  gt_assert(mmsi1 != NULL);
  gt_assert(mmsi2 != NULL);
  return mmsi1->lcpitv.left == mmsi2->lcpitv.left &&
         mmsi1->lcpitv.right == mmsi2->lcpitv.right;
}

void gt_freemmsearchiterator(MMsearchiterator **mmsi)
{
  gt_assert((*mmsi) != NULL);
  gt_encseq_reader_delete((*mmsi)->esr);
  FREESPACE(*mmsi);
}

static bool isleftmaximal(const GtEncseq *dbencseq,
                          GtReadmode readmode,
                          unsigned long dbstart,
                          const Querysubstring *querysubstring)
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
      dbleftchar != accessquery(querysubstring->queryrep,
                                querysubstring->offset-1))
  {
    return true;
  }
  return false;
}

static unsigned long extendright(const GtEncseq *dbencseq,
                          GtEncseqReader *esr,
                          GtReadmode readmode,
                          unsigned long totallength,
                          unsigned long dbend,
                          const Querysubstring *querysubstring,
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
        dbchar != accessquery(querysubstring->queryrep,querypos))
    {
      break;
    }
  }
  return dbpos - dbend;
}

static int runquerysubstringmatch(bool selfmatch,
                                  const GtEncseq *dbencseq,
                                  const ESASuffixptr *suftabpart,
                                  GtReadmode readmode,
                                  unsigned long numberofsuffixes,
                                  uint64_t queryunitnum,
                                  const Queryrep *queryrep,
                                  unsigned long minmatchlength,
                                  Processquerymatch processquerymatch,
                                  void *processquerymatchinfo,
                                  Querymatch *querymatchspaceptr,
                                  GtError *err)
{
  MMsearchiterator *mmsi;
  unsigned long dbstart, totallength, extend;
  uint64_t localqueryunitnum = queryunitnum;
  unsigned long localqueryoffset = 0;
  Querysubstring querysubstring;
  bool haserr = false;

  gt_assert(numberofsuffixes > 0);
  totallength = gt_encseq_total_length(dbencseq);
  querysubstring.queryrep = queryrep;
  for (querysubstring.offset = 0;
       querysubstring.offset <= queryrep->length - minmatchlength;
       querysubstring.offset++)
  {
    mmsi = newmmsearchiterator_generic(dbencseq,
                                       suftabpart,
                                       0, /* leftbound */
                                       numberofsuffixes-1, /* rightbound */
                                       0, /* offset */
                                       readmode,
                                       &querysubstring,
                                       minmatchlength);
    while (!haserr && gt_nextmmsearchiterator(&dbstart,mmsi))
    {
      if (isleftmaximal(dbencseq,
                        readmode,
                        dbstart,
                        &querysubstring))
      {
        extend = extendright(dbencseq,
                             mmsi->esr,
                             readmode,
                             totallength,
                             dbstart + minmatchlength,
                             &querysubstring,
                             minmatchlength);
        gt_querymatch_fill(querymatchspaceptr,
                        extend + minmatchlength,
                        dbstart,
                        queryrep->readmode,
                        selfmatch,
                        localqueryunitnum,
                        localqueryoffset,
                        queryrep->length);
        if (processquerymatch(processquerymatchinfo,
                              dbencseq,
                              querymatchspaceptr,
                              err) != 0)
        {
          haserr = true;
        }
      }
    }
    gt_freemmsearchiterator(&mmsi);
    if (!haserr)
    {
      if (accessquery(queryrep,querysubstring.offset) == (GtUchar) SEPARATOR)
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

int gt_callenumquerymatches(const char *indexname,
                         const GtStrArray *queryfiles,
                         bool echoquery,
                         unsigned int userdefinedleastlength,
                         Processquerymatch processquerymatch,
                         void *processquerymatchinfo,
                         GtLogger *logger,
                         GtError *err)
{
  Suffixarray suffixarray;
  unsigned long totallength = 0;
  bool haserr = false;
  Querymatch *querymatchspaceptr = gt_querymatch_new();

  if (gt_mapsuffixarray(&suffixarray,
                     SARR_ESQTAB | SARR_SUFTAB | SARR_SSPTAB,
                     indexname,
                     logger,
                     err) != 0)
  {
    haserr = true;
  } else
  {
    totallength = gt_encseq_total_length(suffixarray.encseq);
  }
  if (!haserr && echoquery)
  {
    if (gt_echodescriptionandsequence(queryfiles,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    GtSeqIterator *seqit;
    const GtUchar *query;
    unsigned long querylen;
    char *desc = NULL;
    int retval;
    uint64_t queryunitnum;

    seqit = gt_seqiterator_sequence_buffer_new(queryfiles, err);
    if (seqit == NULL)
    {
      haserr = true;
    }
    if (!haserr)
    {
      gt_seqiterator_set_symbolmap(seqit,
                      gt_alphabet_symbolmap(gt_encseq_alphabet(
                                                          suffixarray.encseq)));
      for (queryunitnum = 0; /* Nothing */; queryunitnum++)
      {
        retval = gt_seqiterator_next(seqit,
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
        if (querylen >= (unsigned long) userdefinedleastlength)
        {
          Queryrep queryrep;

          queryrep.sequence = query;
          queryrep.encseq = NULL;
          queryrep.readmode = GT_READMODE_FORWARD;
          queryrep.startpos = 0;
          queryrep.length = (unsigned long) querylen;
          if (runquerysubstringmatch(false,
                                     suffixarray.encseq,
                                     suffixarray.suftab,
                                     suffixarray.readmode,
                                     totallength+1,
                                     queryunitnum,
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
      gt_seqiterator_delete(seqit);
    }
  }
  gt_querymatch_delete(querymatchspaceptr);
  gt_freesuffixarray(&suffixarray);
  return haserr ? -1 : 0;
}

int gt_callenumselfmatches(const char *indexname,
                        GtReadmode queryreadmode,
                        unsigned int userdefinedleastlength,
                        Processquerymatch processquerymatch,
                        void *processquerymatchinfo,
                        GtLogger *logger,
                        GtError *err)
{
  Suffixarray suffixarray;
  unsigned long totallength = 0;
  bool haserr = false;
  Querymatch *querymatchspaceptr = gt_querymatch_new();

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
    totallength = gt_encseq_total_length(suffixarray.encseq);
  }
  if (!haserr)
  {
    unsigned long seqnum, numofsequences, seqlength, seqstartpos;
    Queryrep queryrep;

    numofsequences = gt_encseq_num_of_sequences(suffixarray.encseq);
    queryrep.sequence = NULL;
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
        if (runquerysubstringmatch(true,
                                   suffixarray.encseq,
                                   suffixarray.suftab,
                                   suffixarray.readmode,
                                   totallength+1,
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
  }
  gt_querymatch_delete(querymatchspaceptr);
  gt_freesuffixarray(&suffixarray);
  return haserr ? -1 : 0;
}

static int constructsarrandrunmmsearch(
                 const GtEncseq *dbencseq,
                 GtReadmode readmode,
                 unsigned int prefixlength,
                 unsigned int numofparts,
                 unsigned long maximumspace,
                 const GtUchar *query,
                 unsigned long querylen,
                 unsigned int minlength,
                 Processquerymatch processquerymatch,
                 void *processquerymatchinfo,
                 GtTimer *sfxprogress,
                 bool withprogressbar,
                 GtError *err)
{
  unsigned long numofsuffixes;
  bool haserr = false;
  Sfxiterator *sfi;
  Queryrep queryrep;
  Querymatch *querymatchspaceptr = gt_querymatch_new();
  Sfxstrategy sfxstrategy;

  defaultsfxstrategy(&sfxstrategy,
                     gt_encseq_bitwise_cmp_ok(dbencseq) ? false : true);
  sfxstrategy.outsuftabonfile = false;
  sfi = gt_Sfxiterator_new(dbencseq,
                           readmode,
                           prefixlength,
                           numofparts,
                           maximumspace,
                           NULL, /* outlcpinfo */
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
    queryrep.sequence = query;
    queryrep.encseq = NULL;
    queryrep.readmode = GT_READMODE_FORWARD;
    queryrep.startpos = 0;
    queryrep.length = (unsigned long) querylen;
    while (true)
    {
      suffixsortspace = gt_Sfxiterator_next(&numofsuffixes,NULL,sfi);
      if (suffixsortspace == NULL)
      {
        break;
      }
      if (runquerysubstringmatch(false,
                                 dbencseq,
                                 (const ESASuffixptr *)
                                 gt_suffixsortspace_ulong_get(
                                          suffixsortspace),
                                 readmode,
                                 numofsuffixes,
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
  }
  gt_querymatch_delete(querymatchspaceptr);
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
                               Processquerymatch processquerymatch,
                               void *processquerymatchinfo,
                               GtLogger *logger,
                               GtError *err)
{
  unsigned int numofchars;
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
  if (constructsarrandrunmmsearch(dbencseq,
                                  GT_READMODE_FORWARD,
                                  gt_recommendedprefixlength(numofchars,dblen,
                                                             true),
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
