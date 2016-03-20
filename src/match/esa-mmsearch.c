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
  const GtEncseq *encseq;
  GtReadmode readmode;
  GtUword startpos, seqlen;
} GtQueryrepresentation;

typedef struct
{
  GtQueryrepresentation *queryrep;
  GtUword currentoffset; /* position relative to startpos */
} GtQuerysubstring;

static GtUchar gt_mmsearch_accessquery(const GtQueryrepresentation *queryrep,
                                       GtUword pos)
{
  GtUword abspos, cc;

  gt_assert(queryrep != NULL);
  gt_assert(pos < queryrep->seqlen);
  abspos = queryrep->startpos + (queryrep->readmode == GT_READMODE_FORWARD
                                  ? pos
                                  : GT_REVERSEPOS(queryrep->seqlen,pos));
  if (queryrep->sequence != NULL)
  {
    cc = queryrep->sequence[abspos];
  } else
  {
    gt_assert(queryrep->encseq != NULL);
    cc = gt_encseq_get_encoded_char(queryrep->encseq,abspos,
                                    GT_READMODE_FORWARD);
  }
  if (GT_ISDIRCOMPLEMENT(queryrep->readmode))
  {
    if (ISSPECIAL(cc))
    {
      return cc;
    }
    return GT_COMPLEMENTBASE(cc);
  } else
  {
    return cc;
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
          if ((LCPLEN) >= (GtUword) minmatchlength)\
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
                                         querysubstring->currentoffset + \
                                          (LCPLEN));\
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
                        GtUword minmatchlength)
{
  GtUword left, leftsave, mid, right, lpref, rpref,
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
  GtUword sufindex;
  const ESASuffixptr *suftab;
  GtEncseqReader *esr;
};

static void gt_mmsearchiterator_reinit(GtMMsearchiterator *mmsi,
                                       const GtEncseq *dbencseq,
                                       const ESASuffixptr *suftab,
                                       GtUword leftbound,
                                       GtUword rightbound,
                                       GtUword itvoffset,
                                       GtReadmode readmode,
                                       const GtQuerysubstring *querysubstring,
                                       GtUword minmatchlength)
{
  mmsi->suftab = suftab;
  if (mmsi->esr == NULL)
  {
    mmsi->esr = gt_encseq_create_reader_with_readmode(dbencseq, readmode, 0);
  } else
  {
    gt_encseq_reader_reinit_with_readmode(mmsi->esr,dbencseq,readmode,0);
  }
  mmsi->lcpitv.left = leftbound;
  mmsi->lcpitv.right = rightbound;
  mmsi->lcpitv.offset = itvoffset;
  if (!gt_mmsearch(dbencseq,mmsi->esr,suftab,readmode,&mmsi->lcpitv,
                   querysubstring,minmatchlength))
  {
    mmsi->lcpitv.left = 1UL;
    mmsi->lcpitv.right = 0;
  }
  mmsi->sufindex = mmsi->lcpitv.left;
}

static GtMMsearchiterator *gt_mmsearchiterator_new_empty(void)
{
  GtMMsearchiterator *mmsi = gt_malloc(sizeof *mmsi);

  mmsi->esr = NULL;
  return mmsi;
}

static GtMMsearchiterator *gt_mmsearchiterator_new(
                                       const GtEncseq *dbencseq,
                                       const ESASuffixptr *suftab,
                                       GtUword leftbound,
                                       GtUword rightbound,
                                       GtUword itvoffset,
                                       GtReadmode readmode,
                                       const GtQuerysubstring *querysubstring,
                                       GtUword minmatchlength)
{
  GtMMsearchiterator *mmsi = gt_mmsearchiterator_new_empty();

  gt_mmsearchiterator_reinit(mmsi,
                             dbencseq,
                             suftab,
                             leftbound,
                             rightbound,
                             itvoffset,
                             readmode,
                             querysubstring,
                             minmatchlength);
  return mmsi;
}

GtMMsearchiterator *gt_mmsearchiterator_new_complete_plain(
                                   const GtEncseq *dbencseq,
                                   const void *voidsuftab, /* XXX */
                                   GtUword leftbound,
                                   GtUword rightbound,
                                   GtUword itvoffset,
                                   GtReadmode readmode,
                                   const GtUchar *pattern,
                                   GtUword patternlength)
{
  GtQueryrepresentation queryrep;
  GtQuerysubstring querysubstring;
  const ESASuffixptr *suftab = (const ESASuffixptr *) voidsuftab; /* XXX */

  queryrep.sequence = pattern;
  queryrep.encseq = NULL;
  queryrep.readmode = GT_READMODE_FORWARD;
  queryrep.startpos = 0;
  queryrep.seqlen = patternlength;
  querysubstring.queryrep = &queryrep;
  querysubstring.currentoffset = 0;
  return gt_mmsearchiterator_new(dbencseq,
                                 suftab,
                                 leftbound,
                                 rightbound,
                                 itvoffset,
                                 readmode,
                                 &querysubstring,
                                 patternlength);
}

GtUword gt_mmsearchiterator_count(const GtMMsearchiterator *mmsi)
{
  gt_assert(mmsi != NULL);
  if (mmsi->lcpitv.left > mmsi->lcpitv.right)
  {
    return 0;
  }
  return mmsi->lcpitv.right - mmsi->lcpitv.left + 1;
}

bool gt_mmsearchiterator_next(GtUword *dbstart,GtMMsearchiterator *mmsi)
{
  gt_assert(mmsi != NULL);
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
  gt_assert(mmsi1 != NULL && mmsi2 != NULL);
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
                                      GtUword dbstart,
                                      const GtQuerysubstring *querysubstring)
{
  GtUchar dbleftchar;

  if (dbstart == 0 || querysubstring->currentoffset == 0)
  {
    return true;
  }
  dbleftchar = gt_encseq_get_encoded_char(dbencseq, /* Random access */
                                          dbstart-1,
                                          readmode);
  if (ISSPECIAL(dbleftchar) ||
      dbleftchar != gt_mmsearch_accessquery(querysubstring->queryrep,
                                            querysubstring->currentoffset-1))
  {
    return true;
  }
  return false;
}

static bool gt_mum_isleftmaximal(const GtEncseq *dbencseq,
                                 GtReadmode readmode,
                                 GtUword dbstart,
                                 GtUword queryoffset,
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

static GtUword gt_mmsearch_extendright(const GtEncseq *dbencseq,
                                             GtEncseqReader *esr,
                                             GtReadmode readmode,
                                             GtUword totallength,
                                             GtUword dbend,
                                             const GtQuerysubstring
                                               *querysubstring,
                                             GtUword matchlength)
{
  GtUchar dbchar;
  GtUword dbpos, querypos;

  if (dbend < totallength)
  {
    gt_encseq_reader_reinit_with_readmode(esr,dbencseq,readmode,dbend);
  }
  for (dbpos = dbend, querypos = querysubstring->currentoffset + matchlength;
       dbpos < totallength &&
       querypos < querysubstring->queryrep->seqlen;
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

void gt_queryuniquematch(bool selfmatch,
                        const Suffixarray *suffixarray,
                        uint64_t queryunitnum,
                        GtQueryrepresentation *queryrep,
                        GtUword minmatchlength,
                        GtProcessquerymatch processquerymatch,
                        void *processquerymatchinfo,
                        GtQuerymatch *querymatchspaceptr)
{
  GtUword offset, totallength = gt_encseq_total_length(suffixarray->encseq),
          localqueryoffset = 0;
  uint64_t localqueryunitnum = queryunitnum;

  gt_assert(!selfmatch && queryrep->seqlen >= minmatchlength);
  for (offset = 0; offset <= queryrep->seqlen - minmatchlength; offset++)
  {
    GtUword matchlen, dbstart;

    matchlen = gt_suffixarrayfindmums (suffixarray,
                                       0,
                                       0, /* leftbound */
                                       totallength, /* rightbound */
                                       &dbstart,
                                       queryrep->sequence + offset,
                                       queryrep->sequence + queryrep->seqlen);
    if (dbstart != ULONG_MAX &&
        matchlen >= minmatchlength &&
        gt_mum_isleftmaximal(suffixarray->encseq,
                             suffixarray->readmode,
                             dbstart,
                             offset,
                             queryrep->sequence))
    {
      GtUword dbseqnum = gt_encseq_seqnum(suffixarray->encseq,dbstart),
              dbseqstartpos = gt_encseq_seqstartpos(suffixarray->encseq,
                                                    dbseqnum),
              dbseqlen = gt_encseq_seqlength(suffixarray->encseq,dbseqnum);

      gt_querymatch_init(querymatchspaceptr,
                         matchlen,
                         dbstart,
                         dbseqnum,
                         dbstart - dbseqstartpos,
                         dbseqlen,
                         0, /* score */
                         0, /* edist */
                         selfmatch,
                         localqueryunitnum,
                         matchlen,
                         localqueryoffset,
                         queryrep->seqlen);
      processquerymatch(processquerymatchinfo,querymatchspaceptr);
    }
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

static void gt_querysubstringmatch(bool selfmatch,
                                   const GtEncseq *dbencseq,
                                   const ESASuffixptr *suftabpart,
                                   GtReadmode readmode,
                                   GtUword numberofsuffixes,
                                   uint64_t queryunitnum,
                                   GtQueryrepresentation *queryrep,
                                   GtUword minmatchlength,
                                   GtProcessquerymatch processquerymatch,
                                   void *processquerymatchinfo,
                                   GtQuerymatch *querymatchspaceptr)
{
  GtMMsearchiterator *mmsi;
  GtUword totallength, localqueryoffset = 0;
  uint64_t localqueryunitnum = queryunitnum;
  GtQuerysubstring querysubstring;

  gt_assert(numberofsuffixes > 0);
  totallength = gt_encseq_total_length(dbencseq);
  querysubstring.queryrep = queryrep;
  for (querysubstring.currentoffset = 0;
       querysubstring.currentoffset <= queryrep->seqlen - minmatchlength;
       querysubstring.currentoffset++)
  {
    GtUword dbstart;

    mmsi = gt_mmsearchiterator_new(dbencseq,
                                   suftabpart,
                                   0, /* leftbound */
                                   numberofsuffixes - 1, /* rightbound */
                                   0, /* offset */
                                   readmode,
                                   &querysubstring,
                                   minmatchlength);
    while (gt_mmsearchiterator_next(&dbstart,mmsi))
    {
      if (gt_mmsearch_isleftmaximal(dbencseq,
                                    readmode,
                                    dbstart,
                                    &querysubstring))
      {
        GtUword dbseqnum, dbseqstartpos, dbseqlen, extend;

        extend = gt_mmsearch_extendright(dbencseq,
                                         mmsi->esr,
                                         readmode,
                                         totallength,
                                         dbstart + minmatchlength,
                                         &querysubstring,
                                         minmatchlength);

        if (gt_encseq_has_multiseq_support(dbencseq))
        {
          dbseqnum = gt_encseq_seqnum(dbencseq,dbstart);
          dbseqstartpos = gt_encseq_seqstartpos(dbencseq,dbseqnum);
          dbseqlen = gt_encseq_seqlength(dbencseq,dbseqnum);
        } else
        {
          dbseqnum = dbseqstartpos = dbseqlen = 0;
        }
        gt_querymatch_init(querymatchspaceptr,
                           minmatchlength + extend,
                           dbstart,
                           dbseqnum,
                           dbstart - dbseqstartpos,
                           dbseqlen,
                           0, /* score */
                           0, /* edist */
                           selfmatch,
                           localqueryunitnum,
                           minmatchlength + extend,
                           localqueryoffset,
                           queryrep->seqlen);
        processquerymatch(processquerymatchinfo,querymatchspaceptr);
      }
    }
    gt_mmsearchiterator_delete(mmsi);
    mmsi = NULL;
    if (gt_mmsearch_accessquery(queryrep,querysubstring.currentoffset)
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

static int gt_constructsarrandrunmmsearch(
                 const GtEncseq *dbencseq,
                 GtReadmode readmode,
                 unsigned int prefixlength,
                 unsigned int numofparts,
                 GtUword maximumspace,
                 const GtUchar *query,
                 GtUword querylen,
                 unsigned int minlength,
                 GtProcessquerymatch processquerymatch,
                 void *processquerymatchinfo,
                 GtTimer *sfxprogress,
                 bool withprogressbar,
                 GtLogger *logger,
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
                           logger, /* logger */
                           err);
  if (sfi == NULL)
  {
    haserr = true;
  } else
  {
    const GtSuffixsortspace *suffixsortspace;
    GtUword numberofsuffixes;
    GtQuerymatch *querymatchspaceptr = gt_querymatch_new();
    GtQueryrepresentation queryrep;

    queryrep.sequence = query;
    queryrep.encseq = NULL;
    queryrep.readmode = GT_READMODE_FORWARD;
    queryrep.startpos = 0;
    queryrep.seqlen = querylen;
    while (true)
    {
      suffixsortspace = gt_Sfxiterator_next(&numberofsuffixes,NULL,sfi);
      if (suffixsortspace == NULL)
      {
        break;
      }
      gt_querysubstringmatch(false,
                             dbencseq,
                             (const ESASuffixptr *)
                             gt_suffixsortspace_ulong_get(suffixsortspace),
                             readmode,
                             numberofsuffixes,
                             0,
                             &queryrep,
                             (GtUword) minlength,
                             processquerymatch,
                             processquerymatchinfo,
                             querymatchspaceptr);
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
                               GtUword dblen,
                               const GtUchar *query,
                               GtUword querylen,
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

  gt_assert(querylen >= (GtUword) minlength && dblen >= (GtUword) minlength);
  eb = gt_encseq_builder_new(alpha);
  gt_encseq_builder_disable_multiseq_support(eb);
  gt_encseq_builder_disable_description_support(eb);
  gt_encseq_builder_set_logger(eb, logger);
  gt_encseq_builder_add_multiple_encoded(eb,dbseq,dblen);
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
                                     logger,
                                     err) != 0)
  {
    haserr = true;
  }
  gt_encseq_delete(dbencseq);
  dbencseq = NULL;
  return haserr ? -1 : 0;
}

struct GtQuerysubstringmatchiterator
{
  const GtEncseq *dbencseq;
  const ESASuffixptr *suftabpart;
  GtReadmode db_readmode;
  GtUword numberofsuffixes,
          totallength,
          userdefinedleastlength;
  GtMMsearchiterator *mmsi;
  GtQueryrepresentation queryrep;
  GtQuerysubstring querysubstring;
  const GtUchar *query_for_seqit;
  GtUword query_seqlen;
  GtSeqIterator *seqit;
  GtUword dbstart,
          matchlength;
  uint64_t queryunitnum, query_encseq_numofsequences;
  char *desc;
  bool mmsi_defined;
};

GtQuerysubstringmatchiterator *gt_querysubstringmatchiterator_new(
                                     const GtEncseq *dbencseq,
                                     GtUword totallength,
                                     const ESASuffixptr *suftabpart,
                                     GtReadmode db_readmode,
                                     GtUword numberofsuffixes,
                                     const GtStrArray *query_files,
                                     const GtEncseq *query_encseq,
                                     GtReadmode query_readmode,
                                     unsigned int userdefinedleastlength,
                                     GtError *err)
{
  GtQuerysubstringmatchiterator *qsmi = gt_malloc(sizeof *qsmi);

  qsmi->dbencseq = dbencseq;
  qsmi->suftabpart = suftabpart;
  qsmi->db_readmode = db_readmode;
  qsmi->numberofsuffixes = numberofsuffixes;
  qsmi->totallength = totallength;
  qsmi->userdefinedleastlength = (GtUword) userdefinedleastlength;
  qsmi->queryunitnum = 0;
  qsmi->desc = NULL;
  qsmi->query_for_seqit = NULL;
  qsmi->query_seqlen = 0;
  qsmi->queryrep.sequence = NULL;
  qsmi->queryrep.encseq = query_encseq;
  qsmi->queryrep.readmode = query_readmode;
  qsmi->queryrep.startpos = 0;
  qsmi->dbstart = 0;
  qsmi->matchlength = 0;
  qsmi->querysubstring.queryrep = &qsmi->queryrep;
  qsmi->mmsi = gt_mmsearchiterator_new_empty();
  qsmi->mmsi_defined = false;
  if (query_files == NULL || gt_str_array_size(query_files) == 0)
  {
    gt_assert(query_encseq != NULL);
    qsmi->seqit = NULL;
    qsmi->query_encseq_numofsequences
      = (uint64_t) gt_encseq_num_of_sequences(query_encseq);
  } else
  {
    gt_assert(query_encseq == NULL);
    qsmi->seqit = gt_seq_iterator_sequence_buffer_new(query_files, err);
    if (qsmi->seqit == NULL)
    {
      gt_querysubstringmatchiterator_delete(qsmi);
      return NULL;
    }
    gt_seq_iterator_set_symbolmap(qsmi->seqit,
                        gt_alphabet_symbolmap(gt_encseq_alphabet(dbencseq)));
  }
  return qsmi;
}

void gt_querysubstringmatchiterator_delete(GtQuerysubstringmatchiterator *qsmi)
{
  if (qsmi != NULL)
  {
    gt_mmsearchiterator_delete(qsmi->mmsi);
    gt_seq_iterator_delete(qsmi->seqit);
    gt_free(qsmi);
  }
}

GtUword gt_querysubstringmatchiterator_dbstart(
                      const GtQuerysubstringmatchiterator *qsmi)
{
  gt_assert(qsmi != NULL);
  return qsmi->dbstart;
}

GtUword gt_querysubstringmatchiterator_querystart(
                      const GtQuerysubstringmatchiterator *qsmi)
{
  gt_assert(qsmi != NULL);
  return qsmi->querysubstring.currentoffset;
}

GtUword gt_querysubstringmatchiterator_matchlength(
                      const GtQuerysubstringmatchiterator *qsmi)
{
  gt_assert(qsmi != NULL);
  return qsmi->matchlength;
}

uint64_t gt_querysubstringmatchiterator_queryunitnum(
                      const GtQuerysubstringmatchiterator *qsmi)
{
  gt_assert(qsmi != NULL);
  return qsmi->queryunitnum;
}

GtUword gt_querysubstringmatchiterator_query_seqlen(
                      const GtQuerysubstringmatchiterator *qsmi)
{
  gt_assert(qsmi != NULL);
  return qsmi->query_seqlen;
}

const GtUchar *gt_querysubstringmatchiterator_query(
                      const GtQuerysubstringmatchiterator *qsmi)
{
  gt_assert(qsmi != NULL && qsmi->query_for_seqit != NULL);
  return qsmi->query_for_seqit;
}

int gt_querysubstringmatchiterator_next(GtQuerysubstringmatchiterator *qsmi,
                                        GtError *err)
{
  gt_assert(qsmi != NULL);
  while (true)
  {
    if (qsmi->query_seqlen < qsmi->userdefinedleastlength)
    {
      if (qsmi->seqit != NULL)
      {
        int retval = gt_seq_iterator_next(qsmi->seqit,
                                          &qsmi->query_for_seqit,
                                          &qsmi->query_seqlen,
                                          &qsmi->desc,
                                          err);
        if (retval < 0)
        {
          return -1; /* error */
        }
        if (retval == 0)
        {
          return 1; /* no more sequences */
        }
        gt_assert(qsmi->query_seqlen > 0 && qsmi->query_for_seqit != NULL);
        qsmi->queryrep.sequence = qsmi->query_for_seqit;
      } else
      {
        if (qsmi->queryunitnum == qsmi->query_encseq_numofsequences)
        {
          return 1;
        }
        qsmi->queryrep.startpos = gt_encseq_seqstartpos(qsmi->queryrep.encseq,
                                                        qsmi->queryunitnum);
        qsmi->query_seqlen = gt_encseq_seqlength(qsmi->queryrep.encseq,
                                                 qsmi->queryunitnum);
      }
      gt_assert(qsmi->query_seqlen > 0);
      qsmi->queryrep.seqlen = qsmi->query_seqlen;
      qsmi->querysubstring.currentoffset = 0;
    }
    if (qsmi->query_seqlen >= qsmi->userdefinedleastlength)
    {
      if (!qsmi->mmsi_defined)
      {
        gt_mmsearchiterator_reinit(qsmi->mmsi,
                                   qsmi->dbencseq,
                                   qsmi->suftabpart,
                                   0, /* l */
                                   qsmi->numberofsuffixes - 1, /* r */
                                   0, /* offset */
                                   qsmi->db_readmode,
                                   &qsmi->querysubstring,
                                   qsmi->userdefinedleastlength);
        qsmi->mmsi_defined = true;
      } else
      {
        if (gt_mmsearchiterator_next(&qsmi->dbstart,qsmi->mmsi))
        {
          GtUword extend;

          if (gt_mmsearch_isleftmaximal(qsmi->dbencseq,
                                        qsmi->db_readmode,
                                        qsmi->dbstart,
                                        &qsmi->querysubstring))
          {
            extend = gt_mmsearch_extendright(qsmi->dbencseq,
                                             qsmi->mmsi->esr,
                                             qsmi->db_readmode,
                                             qsmi->totallength,
                                             qsmi->dbstart +
                                               qsmi->userdefinedleastlength,
                                             &qsmi->querysubstring,
                                             qsmi->userdefinedleastlength);
            qsmi->matchlength = qsmi->userdefinedleastlength + extend;
            return 0;
          }
        } else
        {
          qsmi->mmsi_defined = false;
          if (qsmi->querysubstring.currentoffset +
              qsmi->userdefinedleastlength < qsmi->query_seqlen)
          {
            qsmi->querysubstring.currentoffset++;
          } else
          {
            qsmi->query_seqlen = 0;
            qsmi->queryunitnum++;
          }
        }
      }
    } else
    {
      qsmi->query_seqlen = 0;
      qsmi->queryunitnum++;
    }
  }
}
