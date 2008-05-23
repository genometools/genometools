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

#include "libgtcore/chardef.h"
#include "libgtcore/minmax.h"
#include "libgtcore/seqiterator.h"
#include "libgtcore/symboldef.h"
#include "sarr-def.h"
#include "seqpos-def.h"
#include "divmodmul.h"
#include "spacedef.h"
#include "esa-mmsearch-def.h"
#include "sfx-suffixer.h"
#include "measure-time-if.h"
#include "format64.h"
#include "stamp.h"

#include "sfx-apfxlen.pr"
#include "esa-map.pr"
#include "echoseq.pr"

#define COMPARE(OFFSET,LCPLEN)\
        sidx = (OFFSET) + (LCPLEN);\
        if (sidx < totallength)\
        {\
          initEncodedsequencescanstate(esr,dbencseq,readmode,sidx);\
        }\
        for (/* Nothing */ ; /* Nothing */; sidx++, (LCPLEN)++)\
        {\
          if (LCPLEN >= (Seqpos) querylen)\
          {\
            retcode = 0;\
            break;\
          }\
          if (sidx >= totallength)\
          {\
            retcode = -1;\
            break;\
          }\
          currentchar = sequentialgetencodedchar(dbencseq,esr,sidx,readmode);\
          retcode = (int) (query[LCPLEN] - currentchar);\
          if (retcode == 0)\
          {\
            if (ISSPECIAL(currentchar) && ISSPECIAL(query[LCPLEN]))\
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
                     const Uchar *query,
                     unsigned long querylen)
{
  Seqpos left, leftsave, mid, right, lpref, rpref, totallength, lcplen, sidx;
  int retcode = 0;
  Uchar currentchar;

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
        mid = DIV2(left+right);
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
        mid = DIV2(left+right);
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

MMsearchiterator *newmmsearchiterator(const Encodedsequence *dbencseq,
                                      const Seqpos *suftab,
                                      Seqpos leftbound,
                                      Seqpos rightbound,
                                      Seqpos offset,
                                      Readmode readmode,
                                      const Uchar *pattern,
                                      unsigned long patternlen)
{
  MMsearchiterator *mmsi;

  ALLOCASSIGNSPACE(mmsi,NULL,MMsearchiterator,1);
  mmsi->lcpitv.left = leftbound;
  mmsi->lcpitv.right = rightbound;
  mmsi->lcpitv.offset = offset;
  mmsi->suftab = suftab;
  mmsi->esr = newEncodedsequencescanstate();
  if (!mmsearch(dbencseq,mmsi->esr,suftab,readmode,&mmsi->lcpitv,
                pattern,patternlen))
  {
    mmsi->lcpitv.left = (Seqpos) 1;
    mmsi->lcpitv.right = 0;
  }
  mmsi->sufindex = mmsi->lcpitv.left;
  return mmsi;
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
  assert(mmsi1 != NULL);
  assert(mmsi2 != NULL);
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
                          const Uchar *query,
                          unsigned long querystart)
{
  Uchar dbleftchar;

  if (dbstart == 0 || querystart == 0)
  {
    return true;
  }
  dbleftchar = getencodedchar(dbencseq, /* Random access */
                              dbstart-1,
                              readmode);
  if (dbleftchar != query[querystart-1] || ISSPECIAL(dbleftchar))
  {
    return true;
  }
  return false;
}

static unsigned long extendright(const Encodedsequence *dbencseq,
                                 Encodedsequencescanstate *esr,
                                 Readmode readmode,
                                 Seqpos totallength,
                                 Seqpos dbend,
                                 const Uchar *query,
                                 unsigned long queryend,
                                 unsigned long querylength)
{
  Uchar dbchar;
  Seqpos dbpos;
  unsigned long querypos;

  if (dbend < totallength)
  {
    initEncodedsequencescanstate(esr,dbencseq,readmode,dbend);
  }
  for (dbpos = dbend, querypos = queryend; /* Nothing */; dbpos++, querypos++)
  {
    if (dbpos >= totallength || querypos >= querylength)
    {
      break;
    }
    dbchar = sequentialgetencodedchar(dbencseq,esr,dbpos,readmode);
    if (dbchar != query[querypos] || ISSPECIAL(dbchar))
    {
      break;
    }
  }
  return querypos - queryend;
}

int runquerysubstringmatch(const Encodedsequence *dbencseq,
                           const Seqpos *suftabpart,
                           Readmode readmode,
                           Seqpos numberofsuffixes,
                           uint64_t unitnum,
                           const Uchar *query,
                           unsigned long querylen,
                           unsigned int minlength,
                           int (*processmaxmatch)(void *,unsigned long,
                                                  Seqpos,uint64_t,
                                                  unsigned long,Error *),
                           void *processmaxmatchinfo,
                           Error *err)
{
  MMsearchiterator *mmsi;
  Seqpos dbstart, totallength;
  unsigned long extend, currentquerystart;
  uint64_t localunitnum = unitnum;
  unsigned long localqueryoffset = 0;

  assert(numberofsuffixes > 0);
  totallength = getencseqtotallength(dbencseq);
  for (currentquerystart = 0;
       currentquerystart <= querylen - minlength;
       currentquerystart++)
  {
    mmsi = newmmsearchiterator(dbencseq,
                               suftabpart,
                               0, /* leftbound */
                               numberofsuffixes-1, /* rightbound */
                               0, /* offset */
                               readmode,
                               query + currentquerystart,
                               (unsigned long) minlength);
    while (nextmmsearchiterator(&dbstart,mmsi))
    {
      if (isleftmaximal(dbencseq,
                        readmode,
                        dbstart,
                        query,
                        currentquerystart))
      {
        extend = extendright(dbencseq,
                             mmsi->esr,
                             readmode,
                             totallength,
                             dbstart + minlength,
                             query,
                             currentquerystart + minlength,
                             querylen);
        if (processmaxmatch(processmaxmatchinfo,
                            extend + (unsigned long) minlength,
                            dbstart,
                            localunitnum,
                            localqueryoffset,
                            err) != 0)
        {
          return -1;
        }
      }
    }
    freemmsearchiterator(&mmsi);
    if (query[currentquerystart] == (Uchar) SEPARATOR)
    {
      localunitnum++;
      localqueryoffset = 0;
    } else
    {
      localqueryoffset++;
    }
  }
  return 0;
}

int callenumquerymatches(const Str *indexname,
                         const StrArray *queryfiles,
                         bool echoquery,
                         unsigned int userdefinedleastlength,
                         int (*processmaxmatch)(void *,unsigned long,Seqpos,
                                                uint64_t,unsigned long,Error *),
                         void *processmaxmatchinfo,
                         Verboseinfo *verboseinfo,
                         Error *err)
{
  Suffixarray suffixarray;
  Seqpos totallength;
  bool haserr = false;

  if (mapsuffixarray(&suffixarray,
                     &totallength,
                     SARR_ESQTAB | SARR_SUFTAB,
                     indexname,
                     verboseinfo,
                     err) != 0)
  {
    haserr = true;
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
    SeqIterator *seqit;
    const Uchar *query;
    unsigned long querylen;
    char *desc = NULL;
    int retval;
    uint64_t unitnum;

    seqit = seqiterator_new(queryfiles,
                            getsymbolmapAlphabet(suffixarray.alpha),
                            true);
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
      if (runquerysubstringmatch(suffixarray.encseq,
                                 suffixarray.suftab,
                                 suffixarray.readmode,
                                 totallength+1,
                                 unitnum,
                                 query,
                                 querylen,
                                 userdefinedleastlength,
                                 processmaxmatch,
                                 processmaxmatchinfo,
                                 err) != 0)
      {
        haserr = true;
        break;
      }
      FREESPACE(desc);
    }
    seqiterator_delete(seqit);
  }
  freesuffixarray(&suffixarray);
  return haserr ? -1 : 0;
}

static int constructsarrandrunmmsearch(
                 Seqpos specialcharacters,
                 Seqpos realspecialranges,
                 const Encodedsequence *dbencseq,
                 Readmode readmode,
                 unsigned int numofchars,
                 const Uchar *characters,
                 unsigned int prefixlength,
                 unsigned int numofparts,
                 const Uchar *query,
                 unsigned long querylen,
                 unsigned int minlength,
                 int (*processmaxmatch)(void *,unsigned long,Seqpos,
                                        uint64_t,unsigned long,Error *),
                 void *processmaxmatchinfo,
                 Measuretime *mtime,
                 Error *err)
{
  const Seqpos *suftabptr;
  Seqpos numofsuffixes;
  bool haserr = false, specialsuffixes = false;
  Sfxiterator *sfi;

  sfi = newSfxiterator(specialcharacters,
                       realspecialranges,
                       dbencseq,
                       readmode,
                       numofchars,
                       characters,
                       prefixlength,
                       numofparts,
                       NULL, /* outlcpinfo */
                       NULL, /* sfxstrategy */
                       mtime,
                       NULL, /* verboseinfo */
                       err);
  if (sfi == NULL)
  {
    haserr = true;
  } else
  {
    while (true)
    {
      suftabptr = nextSfxiterator(&numofsuffixes,&specialsuffixes,mtime,sfi);
      if (suftabptr == NULL)
      {
        break;
      }
      if (runquerysubstringmatch(dbencseq,
                                suftabptr,
                                readmode,
                                numofsuffixes,
                                0,
                                query,
                                querylen,
                                minlength,
                                processmaxmatch,
                                processmaxmatchinfo,
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

int sarrquerysubstringmatch(const Uchar *dbseq,
                            Seqpos dblen,
                            const Uchar *query,
                            unsigned long querylen,
                            unsigned int minlength,
                            const Alphabet *alpha,
                            int (*processmaxmatch)(void *,unsigned long,Seqpos,
                                                   uint64_t,unsigned long,
                                                   Error *),
                            void *processmaxmatchinfo,
                            Verboseinfo *verboseinfo,
                            Error *err)
{
  Specialcharinfo samplespecialcharinfo;
  unsigned int numofchars;
  bool haserr = false;
  Encodedsequence *dbencseq;

  dbencseq = plain2encodedsequence(true,
                                   &samplespecialcharinfo,
                                   dbseq,
                                   dblen,
                                   NULL,
                                   0,
                                   getmapsizeAlphabet(alpha),
                                   verboseinfo);
  numofchars = getnumofcharsAlphabet(alpha);
  if (constructsarrandrunmmsearch(samplespecialcharinfo.specialcharacters,
                                  samplespecialcharinfo.realspecialranges,
                                  dbencseq,
                                  Forwardmode,
                                  numofchars,
                                  getcharactersAlphabet(alpha),
                                  recommendedprefixlength(numofchars,dblen),
                                  1U, /* parts */
                                  query,
                                  querylen,
                                  minlength,
                                  processmaxmatch,
                                  processmaxmatchinfo,
                                  NULL,
                                  err) != 0)
  {
    haserr = true;
  }
  freeEncodedsequence(&dbencseq);
  return haserr ? -1 : 0;
}
