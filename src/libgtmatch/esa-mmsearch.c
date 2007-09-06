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
#include "stamp.h"

#include "sfx-apfxlen.pr"

#define COMPARE(OFFSET)\
        for (sidx = (OFFSET) + lcplen; /* Nothing */; sidx++, lcplen++)\
        {\
          if (lcplen >= (Seqpos) querylen)\
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
          retcode = (int) (query[lcplen] - currentchar);\
          if (retcode == 0)\
          {\
            if (ISSPECIAL(currentchar) && ISSPECIAL(query[lcplen]))\
            {\
              retcode = (int) -1;\
              break;\
            }\
          } else\
          {\
            break;\
          }\
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
  Seqpos sidx, left, leftsave, mid, right, lpref, rpref, totallength;
  int retcode = 0;
  Uchar currentchar;
  Seqpos lcplen;

  totallength = getencseqtotallength(dbencseq);
  leftsave = left = lcpitv->left;
  right = lcpitv->right;
  lcplen = lcpitv->offset;
  COMPARE(suftab[left]);
  if (retcode > 0)
  {
    lpref = lcplen;
    lcplen = lcpitv->offset;
    COMPARE(suftab[right]);
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
        COMPARE(suftab[mid]);
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
  COMPARE(suftab[left]);
  if (retcode < 0)
  {
    return false;
  } else
  {
    lpref = lcplen;
    lcplen = lcpitv->offset;
    COMPARE(suftab[right]);
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
        COMPARE(suftab[mid]);
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
  mmsi->lcpitv.offset = offset;
  mmsi->lcpitv.left = leftbound;
  mmsi->lcpitv.right = rightbound;
  mmsi->suftab = suftab;
  if (!mmsearch(dbencseq,suftab,readmode,&mmsi->lcpitv,pattern,patternlen))
  {
    mmsi->lcpitv.left = (Seqpos) 1;
    mmsi->lcpitv.right = 0;
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

static int processsuftab(const Encodedsequence *dbencseq,
                         const Seqpos *suftabpart,
                         Readmode readmode,
                         Seqpos numberofsuffixes,
                         const Uchar *query,
                         unsigned long querylen,
                         unsigned int minlength,
                         int (*processmaxmatch)(void *,unsigned long,
                                                Seqpos,unsigned long),
                         void *processmaxmatchinfo,
                         Env *env)
{
  MMsearchiterator *mmsi;
  Seqpos dbstart, totallength;
  unsigned long extend, currentquerystart;

  assert(numberofsuffixes > 0);
  totallength = getencseqtotallength(dbencseq);
  for (currentquerystart = 0;
       currentquerystart <= querylen - minlength;
       currentquerystart++)
  {
    mmsi = newmmsearchiterator(dbencseq,
                               suftabpart,
                               0,
                               0,
                               numberofsuffixes-1,
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
                            currentquerystart) != 0)
        {
          return -1;
	}
      }
    }
    freemmsearchiterator(&mmsi,env);
  }
  return 0;
}

static int suffixeratormatch(
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
                                        unsigned long),
                 void *processmaxmatchinfo,
                 Measuretime *mtime,
                 Env *env)
{
  const Seqpos *suftabptr;
  Seqpos len;
  bool haserr = false;
  Sfxiterator *sfi;

  sfi = newsfxiterator(specialcharacters,
                       specialranges,
                       dbencseq,
                       readmode,
                       numofchars,
                       prefixlength,
                       numofparts,
                       mtime,
                       env);
  if(sfi == NULL)
  {
    haserr = true;
  } else
  {
    while(true)
    {
      suftabptr = nextSfxiterator(&len,mtime,sfi,env);
      if(suftabptr == NULL)
      {
        break;
      }
      if(processsuftab(dbencseq,suftabptr,readmode,len,query,querylen,
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
  if(sfi != NULL)
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
                                                   unsigned long),
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
                                   alpha,
                                   env);
  numofchars = getnumofcharsAlphabet(alpha);
  if (suffixeratormatch(samplespecialcharinfo.specialcharacters,
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
