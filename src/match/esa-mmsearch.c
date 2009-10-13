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
#include "core/seqiterator.h"
#include "core/symboldef.h"
#include "sarr-def.h"
#include "seqpos-def.h"
#include "spacedef.h"
#include "esa-mmsearch.h"
#include "sfx-suffixer.h"
#include "sfx-progress.h"
#include "format64.h"
#include "stamp.h"
#include "esa-map.h"
#include "echoseq.h"
#include "sfx-apfxlen.h"

typedef struct
{
  const GtUchar *sequence;
  const Encodedsequence *encseq;
  Readmode readmode;
  Seqpos startpos, length;
} Queryrep;

typedef struct
{
  const Queryrep *queryrep;
  Seqpos offset; /* position relative to startpos */
} Querysubstring;

static GtUchar accessquery(int line,const Queryrep *queryrep,
                           Seqpos pos)
{
  Seqpos abspos = queryrep->startpos + pos;

  if (pos >= queryrep->length)
  {
    fprintf(stderr,"line %d: pos = " FormatSeqpos " >= "
                   FormatSeqpos "=queryrep->length\n",
                   line,PRINTSeqposcast(pos),
                   PRINTSeqposcast(queryrep->length));
    exit(EXIT_FAILURE);
  }
  gt_assert(pos < queryrep->length);
  if (queryrep->sequence != NULL)
  {
    gt_assert(queryrep->readmode == Forwardmode);
    return queryrep->sequence[abspos];
  } else
  {
    gt_assert(queryrep->readmode != Forwardmode);
    gt_assert(queryrep->encseq != NULL);
    return getencodedchar(queryrep->encseq,abspos,queryrep->readmode);
  }
}

#define COMPARE(OFFSET,LCPLEN)\
        sidx = (OFFSET) + (LCPLEN);\
        if (sidx < totallength)\
        {\
          initEncodedsequencescanstate(esr,dbencseq,readmode,sidx);\
        }\
        for (/* Nothing */ ; /* Nothing */; sidx++, (LCPLEN)++)\
        {\
          if ((LCPLEN) >= (Seqpos) minmatchlength)\
          {\
            retcode = 0;\
            break;\
          }\
          if (sidx >= totallength)\
          {\
            retcode = -1;\
            break;\
          }\
          currentdbchar = sequentialgetencodedchar(dbencseq,esr,sidx,readmode);\
          currentquerychar = accessquery(__LINE__,querysubstring->queryrep,\
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

static bool mmsearch(const Encodedsequence *dbencseq,
                     Encodedsequencescanstate *esr,
                     const Seqpos *suftab,
                     Readmode readmode,
                     Lcpinterval *lcpitv,
                     const Querysubstring *querysubstring,
                     unsigned long minmatchlength)
{
  Seqpos left, leftsave, mid, right, lpref, rpref, totallength, lcplen, sidx;
  int retcode = 0;
  GtUchar currentdbchar, currentquerychar;

  totallength = getencseqtotallength(dbencseq);
  leftsave = left = lcpitv->left;
  right = lcpitv->right;
  lcplen = lcpitv->offset;
  COMPARE(suftab[left],lcplen);
  if (retcode > 0)
  {
    lpref = lcplen;
    lcplen = lcpitv->offset;
    COMPARE(suftab[right],lcplen);
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
        COMPARE(suftab[mid],lcplen);
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
  COMPARE(suftab[left],lcplen);
  if (retcode < 0)
  {
    return false;
  } else
  {
    lpref = lcplen;
    lcplen = lcpitv->offset;
    COMPARE(suftab[right],lcplen);
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
        COMPARE(suftab[mid],lcplen);
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
  Seqpos sufindex;
  const Seqpos *suftab;
  Encodedsequencescanstate *esr;
};

static MMsearchiterator *newmmsearchiterator_generic(
                                       const Encodedsequence *dbencseq,
                                       const Seqpos *suftab,
                                       Seqpos leftbound,
                                       Seqpos rightbound,
                                       Seqpos itvoffset,
                                       Readmode readmode,
                                       const Querysubstring *querysubstring,
                                       unsigned long minmatchlength)

{
  MMsearchiterator *mmsi;

  ALLOCASSIGNSPACE(mmsi,NULL,MMsearchiterator,1);
  mmsi->lcpitv.left = leftbound;
  mmsi->lcpitv.right = rightbound;
  mmsi->lcpitv.offset = itvoffset;
  mmsi->suftab = suftab;
  mmsi->esr = newEncodedsequencescanstate();
  if (!mmsearch(dbencseq,mmsi->esr,suftab,readmode,&mmsi->lcpitv,
                querysubstring,minmatchlength))
  {
    mmsi->lcpitv.left = (Seqpos) 1;
    mmsi->lcpitv.right = 0;
  }
  mmsi->sufindex = mmsi->lcpitv.left;
  return mmsi;
}

MMsearchiterator *newmmsearchiteratorcomplete_plain(
                                   const Encodedsequence *dbencseq,
                                   const Seqpos *suftab,
                                   Seqpos leftbound,
                                   Seqpos rightbound,
                                   Seqpos itvoffset,
                                   Readmode readmode,
                                   const GtUchar *pattern,
                                   unsigned long patternlength)
{
  Queryrep queryrep;
  Querysubstring querysubstring;

  queryrep.sequence = pattern;
  queryrep.encseq = NULL;
  queryrep.readmode = Forwardmode;
  queryrep.startpos = 0;
  queryrep.length = (Seqpos) patternlength;
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

Seqpos countmmsearchiterator(const MMsearchiterator *mmsi)
{
  if (mmsi->lcpitv.left > mmsi->lcpitv.right)
  {
    return 0;
  }
  return mmsi->lcpitv.right - mmsi->lcpitv.left + 1;
}

bool nextmmsearchiterator(Seqpos *dbstart,MMsearchiterator *mmsi)
{
  if (mmsi->sufindex <= mmsi->lcpitv.right)
  {
    *dbstart = mmsi->suftab[mmsi->sufindex++];
    return true;
  }
  return false;
}

bool isemptymmsearchiterator(const MMsearchiterator *mmsi)
{
  return mmsi == NULL || mmsi->lcpitv.left > mmsi->lcpitv.right;
}

bool identicalmmsearchiterators(const MMsearchiterator *mmsi1,
                                const MMsearchiterator *mmsi2)
{
  gt_assert(mmsi1 != NULL);
  gt_assert(mmsi2 != NULL);
  return mmsi1->lcpitv.left == mmsi2->lcpitv.left &&
         mmsi1->lcpitv.right == mmsi2->lcpitv.right;
}

void freemmsearchiterator(MMsearchiterator **mmsi)
{
  freeEncodedsequencescanstate(&(*mmsi)->esr);
  FREESPACE(*mmsi);
}

static bool isleftmaximal(const Encodedsequence *dbencseq,
                          Readmode readmode,
                          Seqpos dbstart,
                          const Querysubstring *querysubstring)
{
  GtUchar dbleftchar;

  if (dbstart == 0 || querysubstring->offset == 0)
  {
    return true;
  }
  dbleftchar = getencodedchar(dbencseq, /* Random access */
                              dbstart-1,
                              readmode);
  if (ISSPECIAL(dbleftchar) ||
      dbleftchar != accessquery(__LINE__,querysubstring->queryrep,
                                querysubstring->offset-1))
  {
    return true;
  }
  return false;
}

static Seqpos extendright(const Encodedsequence *dbencseq,
                          Encodedsequencescanstate *esr,
                          Readmode readmode,
                          Seqpos totallength,
                          Seqpos dbend,
                          const Querysubstring *querysubstring,
                          unsigned long matchlength)
{
  GtUchar dbchar;
  Seqpos dbpos, querypos;

  if (dbend < totallength)
  {
    initEncodedsequencescanstate(esr,dbencseq,readmode,dbend);
  }
  for (dbpos = dbend, querypos = querysubstring->offset + matchlength;
       dbpos < totallength &&
       querypos < querysubstring->queryrep->length;
       dbpos++, querypos++)
  {
    dbchar = sequentialgetencodedchar(dbencseq,esr,dbpos,readmode);
    if (ISSPECIAL(dbchar) ||
        dbchar != accessquery(__LINE__,querysubstring->queryrep,querypos))
    {
      break;
    }
  }
  return dbpos - dbend;
}

static int runquerysubstringmatch(const Encodedsequence *dbencseq,
                                  const Seqpos *suftabpart,
                                  Readmode readmode,
                                  Seqpos numberofsuffixes,
                                  uint64_t queryunitnum,
                                  const Queryrep *queryrep,
                                  unsigned long minmatchlength,
                                  Processquerymatch processquerymatch,
                                  void *processquerymatchinfo,
                                  GtError *err)
{
  MMsearchiterator *mmsi;
  Seqpos dbstart, totallength, extend;
  uint64_t localqueryunitnum = queryunitnum;
  Seqpos localqueryoffset = 0;
  Querysubstring querysubstring;

  gt_assert(numberofsuffixes > 0);
  totallength = getencseqtotallength(dbencseq);
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
    while (nextmmsearchiterator(&dbstart,mmsi))
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
        if (processquerymatch(processquerymatchinfo,
                              dbencseq,
                              extend + minmatchlength,
                              dbstart,
                              queryrep->readmode,
                              localqueryunitnum,
                              localqueryoffset,
                              queryrep->length,
                              err) != 0)
        {
          return -1;
        }
      }
    }
    freemmsearchiterator(&mmsi);
    if (accessquery(__LINE__,queryrep,querysubstring.offset)
        == (GtUchar) SEPARATOR)
    {
      localqueryunitnum++;
      localqueryoffset = 0;
    } else
    {
      localqueryoffset++;
    }
  }
  return 0;
}

int callenumquerymatches(const GtStr *indexname,
                         const GtStrArray *queryfiles,
                         bool echoquery,
                         unsigned int userdefinedleastlength,
                         Processquerymatch processquerymatch,
                         void *processquerymatchinfo,
                         Verboseinfo *verboseinfo,
                         GtError *err)
{
  Suffixarray suffixarray;
  Seqpos totallength = 0;
  bool haserr = false;

  if (mapsuffixarray(&suffixarray,
                     SARR_ESQTAB | SARR_SUFTAB | SARR_SSPTAB,
                     indexname,
                     verboseinfo,
                     err) != 0)
  {
    haserr = true;
  } else
  {
    totallength = getencseqtotallength(suffixarray.encseq);
  }
  if (!haserr && echoquery)
  {
    if (echodescriptionandsequence(queryfiles,err) != 0)
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

    seqit = gt_seqiterator_new(queryfiles, err);
    if (seqit == NULL)
    {
      haserr = true;
    }
    if (!haserr)
    {
      gt_seqiterator_set_symbolmap(seqit,
                                   getencseqAlphabetsymbolmap(suffixarray.
                                                              encseq));
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
          queryrep.readmode = Forwardmode;
          queryrep.startpos = 0;
          queryrep.length = (Seqpos) querylen;
          if (runquerysubstringmatch(suffixarray.encseq,
                                     suffixarray.suftab,
                                     suffixarray.readmode,
                                     totallength+1,
                                     queryunitnum,
                                     &queryrep,
                                     (unsigned long) userdefinedleastlength,
                                     processquerymatch,
                                     processquerymatchinfo,
                                     err) != 0)
          {
            haserr = true;
            break;
          }
        }
        gt_free(desc);
      }
      gt_seqiterator_delete(seqit);
    }
  }
  freesuffixarray(&suffixarray);
  return haserr ? -1 : 0;
}

int callenumselfmatches(const GtStr *indexname,
                        Readmode queryreadmode,
                        unsigned int userdefinedleastlength,
                        Processquerymatch processquerymatch,
                        void *processquerymatchinfo,
                        Verboseinfo *verboseinfo,
                        GtError *err)
{
  Suffixarray suffixarray;
  Seqpos totallength = 0;
  bool haserr = false;

  gt_assert(queryreadmode != Forwardmode);
  if (mapsuffixarray(&suffixarray,
                     SARR_ESQTAB | SARR_SUFTAB | SARR_SSPTAB,
                     indexname,
                     verboseinfo,
                     err) != 0)
  {
    haserr = true;
  } else
  {
    totallength = getencseqtotallength(suffixarray.encseq);
  }
  if (!haserr)
  {
    unsigned long seqnum, numofsequences;
    Queryrep queryrep;
    Seqinfo seqinfo;

    numofsequences = getencseqnumofdbsequences(suffixarray.encseq);
    queryrep.sequence = NULL;
    queryrep.encseq = suffixarray.encseq;
    queryrep.readmode = queryreadmode;
    for (seqnum = 0; seqnum < numofsequences; seqnum++)
    {
      getencseqSeqinfo(&seqinfo,suffixarray.encseq,seqnum);
      if (seqinfo.seqlength >= (Seqpos) userdefinedleastlength)
      {
        queryrep.startpos = seqinfo.seqstartpos;
        queryrep.length = seqinfo.seqlength;
        if (runquerysubstringmatch(suffixarray.encseq,
                                   suffixarray.suftab,
                                   suffixarray.readmode,
                                   totallength+1,
                                   (uint64_t) seqnum,
                                   &queryrep,
                                   (unsigned long) userdefinedleastlength,
                                   processquerymatch,
                                   processquerymatchinfo,
                                   err) != 0)
        {
          haserr = true;
          break;
        }
      }
    }
  }
  freesuffixarray(&suffixarray);
  return haserr ? -1 : 0;
}

static int constructsarrandrunmmsearch(
                 const Encodedsequence *dbencseq,
                 Readmode readmode,
                 unsigned int prefixlength,
                 unsigned int numofparts,
                 const GtUchar *query,
                 unsigned long querylen,
                 unsigned int minlength,
                 Processquerymatch processquerymatch,
                 void *processquerymatchinfo,
                 Sfxprogress *sfxprogress,
                 GtError *err)
{
  const Seqpos *suftabptr;
  Seqpos numofsuffixes;
  bool haserr = false, specialsuffixes = false;
  Sfxiterator *sfi;
  Queryrep queryrep;

  sfi = newSfxiterator(dbencseq,
                       readmode,
                       prefixlength,
                       numofparts,
                       NULL, /* outlcpinfo */
                       NULL, /* sfxstrategy */
                       sfxprogress,
                       NULL, /* verboseinfo */
                       err);
  if (sfi == NULL)
  {
    haserr = true;
  } else
  {
    queryrep.sequence = query;
    queryrep.encseq = NULL;
    queryrep.readmode = Forwardmode;
    queryrep.startpos = 0;
    queryrep.length = (Seqpos) querylen;
    while (true)
    {
      suftabptr = nextSfxiterator(&numofsuffixes,&specialsuffixes,sfi);
      if (suftabptr == NULL)
      {
        break;
      }
      if (runquerysubstringmatch(dbencseq,
                                suftabptr,
                                readmode,
                                numofsuffixes,
                                0,
                                &queryrep,
                                (unsigned long) minlength,
                                processquerymatch,
                                processquerymatchinfo,
                                err) != 0)
      {
        haserr = true;
        break;
      }
    }
  }
  if (sfi != NULL)
  {
    freeSfxiterator(&sfi);
  }
  return haserr ? -1 : 0;
}

int sarrquerysubstringmatch(const GtUchar *dbseq,
                            Seqpos dblen,
                            const GtUchar *query,
                            unsigned long querylen,
                            unsigned int minlength,
                            const GtAlphabet *alpha,
                            Processquerymatch processquerymatch,
                            void *processquerymatchinfo,
                            Verboseinfo *verboseinfo,
                            GtError *err)
{
  unsigned int numofchars;
  bool haserr = false;
  Encodedsequence *dbencseq;

  dbencseq = plain2encodedsequence(true,
                                   dbseq,
                                   dblen,
                                   NULL,
                                   0,
                                   alpha,
                                   verboseinfo);
  numofchars = gt_alphabet_num_of_chars(alpha);
  if (constructsarrandrunmmsearch(dbencseq,
                                  Forwardmode,
                                  recommendedprefixlength(numofchars,dblen),
                                  1U, /* parts */
                                  query,
                                  querylen,
                                  minlength,
                                  processquerymatch,
                                  processquerymatchinfo,
                                  NULL,
                                  err) != 0)
  {
    haserr = true;
  }
  removealpharef(dbencseq);
  encodedsequence_free(&dbencseq);
  return haserr ? -1 : 0;
}
