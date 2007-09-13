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

#include "libgtcore/minmax.h"
#include "sarr-def.h"
#include "seqpos-def.h"
#include "divmodmul.h"
#include "symboldef.h"
#include "spacedef.h"
#include "chardef.h"
#include "esa-mmsearch-def.h"
#include "sfx-suffixer.h"
#include "measure-time-if.h"
#include "iterseq.h"
#include "format64.h"
#include "stamp.h"

#include "sfx-apfxlen.pr"
#include "sfx-map.pr"

/*
#define COMPARE(OFFSET,LCPLEN)\
        for (sidx = (OFFSET) + LCPLEN; ; sidx++, LCPLEN++)\
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
          currentchar = getencodedchar(dbencseq,sidx,readmode);\
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
*/

#define COMPARE(OFFSET,LCPLEN)\
        retcode = comparecharacters(OFFSET,&lcplen,dbencseq,readmode,\
                                    totallength,query,querylen)

static int comparecharacters(Seqpos start,
                             Seqpos *lcplen,
                             const Encodedsequence *dbencseq,
                             Readmode readmode,
                             Seqpos totallength,
                             const Uchar *query,
                             unsigned long querylen)
{
  Seqpos sidx;
  int retcode = 0;
  Uchar currentchar;

  for (sidx = start + *lcplen; /* Nothing */; sidx++, (*lcplen)++)
  {
    if (*lcplen >= (Seqpos) querylen)
    {
      retcode = 0;
      break;
    }
    if (sidx >= totallength)
    {
      retcode = -1;
      break;
    }
    /* printf("getencodedchar at %u\n",sidx); */
    currentchar = getencodedchar(dbencseq,sidx,readmode);
    retcode = (int) (query[*lcplen] - currentchar);
    if (retcode == 0)
    {
      if (ISSPECIAL(currentchar) && ISSPECIAL(query[*lcplen]))
      {
        retcode = -1;
        break;
      }
    } else
    {
      break;
    }
  }
  return retcode;
}

typedef struct
{
  Seqpos offset,
         left,
         right;
} Lcpinterval;

static bool mmsearch(const Encodedsequence *dbencseq,
                     const Seqpos *suftab,
                     Readmode readmode,
                     Lcpinterval *lcpitv,
                     const Uchar *query,
                     unsigned long querylen)
{
  Seqpos left, leftsave, mid, right, lpref, rpref, totallength;
  int retcode = 0;
  Seqpos lcplen;

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
        /* printf("mid=%u,lcplen=%u\n",mid,lcplen); */
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
};

MMsearchiterator *newmmsearchiterator(const Encodedsequence *dbencseq,
                                      const Seqpos *suftab,
                                      Seqpos leftbound,
                                      Seqpos rightbound,
                                      Seqpos offset,
                                      Readmode readmode,
                                      const Uchar *pattern,
                                      unsigned long patternlen,
                                      Env *env)
{
  MMsearchiterator *mmsi;

  ALLOCASSIGNSPACE(mmsi,NULL,MMsearchiterator,1);
  mmsi->lcpitv.left = leftbound;
  mmsi->lcpitv.right = rightbound;
  mmsi->lcpitv.offset = offset;
  mmsi->suftab = suftab;
  if (!mmsearch(dbencseq,suftab,readmode,&mmsi->lcpitv,pattern,patternlen))
  {
    mmsi->lcpitv.left = (Seqpos) 1;
    mmsi->lcpitv.right = 0;
  } else
  {
    /*
    printf("matching interval %u %u\n",mmsi->lcpitv.left,mmsi->lcpitv.right);
    */
  }
  mmsi->sufindex = mmsi->lcpitv.left;
  return mmsi;
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

void freemmsearchiterator(MMsearchiterator **mmsi,Env *env)
{
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
  dbleftchar = getencodedchar(dbencseq,dbstart-1,readmode);
  if (dbleftchar != query[querystart-1] || ISSPECIAL(dbleftchar))
  {
    return true;
  }
  return false;
}

static unsigned long extendright(const Encodedsequence *dbencseq,
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

  for (dbpos = dbend, querypos = queryend; /* Nothing */; dbpos++, querypos++)
  {
    if (dbpos >= totallength || querypos >= querylength)
    {
      break;
    }
    dbchar = getencodedchar(dbencseq,dbpos,readmode);
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
                                                  unsigned long,Env *),
                           void *processmaxmatchinfo,
                           Env *env)
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
                               (unsigned long) minlength,
                               env);
    while (nextmmsearchiterator(&dbstart,mmsi))
    {
      if (isleftmaximal(dbencseq,
                        readmode,
                        dbstart,
                        query,
                        currentquerystart))
      {
        extend = extendright(dbencseq,
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
                            env) != 0)
        {
          return -1;
	}
      }
    }
    freemmsearchiterator(&mmsi,env);
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

static int echothesequence(const StrArray *queryfiles,Env *env)
{
  Scansequenceiterator *sseqit;
  char *desc = NULL;
  const Uchar *sequence;
  unsigned long seqlen;
  bool haserr = false;
  int retval;

  sseqit = newScansequenceiterator(queryfiles,NULL,true,env);
  while (true)
  {
    retval = nextScansequenceiterator(&sequence,
                                      &seqlen,
                                      &desc,
                                      sseqit,
                                      env);
    if (retval < 0)
    {
      haserr = true;
      break;
    }
    if (retval == 0)
    {
      break;
    }
    fastasymbolstringgeneric(stdout,desc,NULL,sequence,seqlen,
                             (unsigned long) 70);
    FREESPACE(desc);
  }
  freeScansequenceiterator(&sseqit,env);
  return haserr ? -1 : 0;
}

int callenumquerymatches(const Str *indexname,
                         const StrArray *queryfiles,
                         bool echoquery,
                         unsigned int userdefinedleastlength,
                         int (*processmaxmatch)(void *,unsigned long,Seqpos,
                                                uint64_t,unsigned long,Env *),
                         void *processmaxmatchinfo,
                         Env *env)
{
  Suffixarray suffixarray;
  Seqpos totallength;
  bool haserr = false;

  if (mapsuffixarray(&suffixarray,
                     &totallength,
                     SARR_ESQTAB | SARR_SUFTAB,
                     indexname,
                     false,
                     env) != 0)
  {
    haserr = true;
  }
  if (!haserr && echoquery)
  {
    if (echothesequence(queryfiles,env) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    Scansequenceiterator *sseqit;
    const Uchar *query;
    unsigned long querylen;
    char *desc = NULL;
    int retval;
    uint64_t unitnum;

    sseqit = newScansequenceiterator(queryfiles,
                                     getsymbolmapAlphabet(suffixarray.alpha),
                                     true,
                                     env);
    for (unitnum = 0; /* Nothing */; unitnum++)
    {
      retval = nextScansequenceiterator(&query,
                                        &querylen,
                                        &desc,
                                        sseqit,
                                        env);
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
                                env) != 0)
      {
        haserr = true;
        break;
      }
      FREESPACE(desc);
    }
    freeScansequenceiterator(&sseqit,env);
  }
  freesuffixarray(&suffixarray,env);
  return haserr ? -1 : 0;
}

static int constructsarrandrunmmsearch(
                 Seqpos specialcharacters,
                 Seqpos specialranges,
                 const Encodedsequence *dbencseq,
                 Readmode readmode,
                 unsigned int numofchars,
                 unsigned int prefixlength,
                 unsigned int numofparts,
                 const Uchar *query,
                 unsigned long querylen,
                 unsigned int minlength,
                 int (*processmaxmatch)(void *,unsigned long,Seqpos,
                                        uint64_t,unsigned long,Env *),
                 void *processmaxmatchinfo,
                 Measuretime *mtime,
                 Env *env)
{
  const Seqpos *suftabptr;
  Seqpos numofsuffixes;
  bool haserr = false, specialsuffixes = false;
  Sfxiterator *sfi;

  sfi = newSfxiterator(specialcharacters,
                       specialranges,
                       dbencseq,
                       readmode,
                       numofchars,
                       prefixlength,
                       numofparts,
                       mtime,
                       env);
  if (sfi == NULL)
  {
    haserr = true;
  } else
  {
    while (true)
    {
      suftabptr = nextSfxiterator(&numofsuffixes,&specialsuffixes,
                                  mtime,sfi,env);
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
                                env) != 0)
      {
        haserr = true;
        break;
      }
    }
  }
  if (sfi != NULL)
  {
    freeSfxiterator(&sfi,env);
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
                                                   Env *),
                            void *processmaxmatchinfo,
                            Env *env)
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
                                   env);
  numofchars = getnumofcharsAlphabet(alpha);
  if (constructsarrandrunmmsearch(samplespecialcharinfo.specialcharacters,
                                  samplespecialcharinfo.specialranges,
                                  dbencseq,
                                  Forwardmode,
                                  numofchars,
                                  recommendedprefixlength(numofchars,dblen),
                                  (unsigned int) 1, /* parts */
                                  query,
                                  querylen,
                                  minlength,
                                  processmaxmatch,
                                  processmaxmatchinfo,
                                  NULL,
		                  env) != 0)
  {
    haserr = true;
  }
  freeEncodedsequence(&dbencseq,env);
  return haserr ? -1 : 0;
}
