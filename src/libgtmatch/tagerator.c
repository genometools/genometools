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

#include <limits.h>
#include "libgtcore/unused.h"
#include "libgtcore/strarray.h"
#include "libgtcore/ma.h"
#include "libgtcore/error.h"
#include "libgtcore/fileutils.h"
#include "libgtcore/seqiterator.h"
#include "libgtcore/arraydef.h"
#include "tagerator.h"
#include "sarr-def.h"
#include "esa-mmsearch-def.h"
#include "intbits.h"
#include "alphadef.h"
#include "esa-myersapm.h"
#include "esa-limdfs.h"
#include "format64.h"

#include "echoseq.pr"
#include "esa-map.pr"

#define MAXTAGSIZE INTWORDSIZE

static void exactpatternmatching(const Encodedsequence *encseq,
                                 Readmode readmode,
                                 const Seqpos *suftab,
                                 Seqpos totallength,
                                 const Uchar *pattern,
                                 unsigned long patternlength,
                                 void (*processmatch)(void *,Seqpos,Seqpos),
                                 void *processmatchinfo)
{
  MMsearchiterator *mmsi;
  Seqpos dbstartpos;

  mmsi = newmmsearchiterator(encseq,
                             suftab,
                             0,  /* leftbound */
                             totallength, /* rightbound */
                             0, /* offset */
                             readmode,
                             pattern,
                             patternlength);
  while (nextmmsearchiterator(&dbstartpos,mmsi))
  {
    processmatch(processmatchinfo,dbstartpos,(Seqpos) patternlength);
  }
  freemmsearchiterator(&mmsi);
}

static void showmatch(UNUSED void *processinfo,Seqpos startpos,Seqpos len)
{
  printf("match " FormatSeqpos " " FormatSeqpos "\n",
          PRINTSeqposcast(startpos),
          PRINTSeqposcast(len));
}

DECLAREARRAYSTRUCT(Seqpos);

static void storematch(void *processinfo,Seqpos startpos,UNUSED Seqpos len)
{
  ArraySeqpos *storetab = (ArraySeqpos *) processinfo;

  STOREINARRAY(storetab,Seqpos,32,startpos);
}

static int cmpdescend(const void *a,const void *b)
{
  Seqpos *valuea = (Seqpos *) a;
  Seqpos *valueb = (Seqpos *) b;

  if (*valuea < *valueb)
  {
    return 1;
  }
  if (*valuea > *valueb)
  {
    return -1;
  }
  return 0;
}

static int dotransformtag(Uchar *transformedtag,
                          const Uchar *symbolmap,
                          const Uchar *currenttag,
                          unsigned long taglen,
                          uint64_t tagnumber,
                          bool replacewildcard,
                          Error *err)
{
  unsigned long idx;
  Uchar charcode;

  if (taglen > (unsigned long) MAXTAGSIZE)
  {
    error_set(err,"tag of length %lu; tags must not be longer than %lu",
                   taglen, (unsigned long) MAXTAGSIZE);
    return -1;
  }
  for (idx = 0; idx < taglen; idx++)
  {
    charcode = symbolmap[currenttag[idx]];
    if (charcode == (Uchar) UNDEFCHAR)
    {
      error_set(err,"undefined character '%c' in tag number " Formatuint64_t,
                currenttag[idx],
                PRINTuint64_tcast(tagnumber));
      return -1;
    }
    if (charcode == (Uchar) WILDCARD)
    {
      if (replacewildcard)
      {
        charcode = 0; /* (Uchar) (drand48() * (mapsize-1)); */
      } else
      {
        error_set(err,"wildcard in tag number " Formatuint64_t,
                  PRINTuint64_tcast(tagnumber));
        return -1;
      }
    }
    transformedtag[idx] = charcode;
  }
  return 0;
}

static void performthesearch(const TageratorOptions *tageratoroptions,
                             Myersonlineresources *mor,
                             Limdfsresources *limdfsresources,
                             const Encodedsequence *encseq,
                             Readmode readmode,
                             const Seqpos *suftab,
                             Seqpos totallength,
                             const Uchar *transformedtag,
                             unsigned long taglen,
                             void (*processmatch)(void *,Seqpos,Seqpos),
                             void *processmatchinfooffline,
                             UNUSED bool rcmatch)
{
  if (tageratoroptions->online || tageratoroptions->docompare)
  {
    edistmyersbitvectorAPM(mor,
                           transformedtag,
                           taglen,
                           tageratoroptions->maxdistance);
  }
  if (!tageratoroptions->online || tageratoroptions->docompare)
  {
    if (tageratoroptions->maxdistance == 0)
    {
      exactpatternmatching(encseq,
                           readmode,
                           suftab,
                           totallength,
                           transformedtag,
                           taglen,
                           processmatch,
                           processmatchinfooffline);
    } else
    {
      esalimiteddfs(limdfsresources,
                    transformedtag,
                    taglen,
                    tageratoroptions->maxdistance);
    }
  }
}

static void compareresults(UNUSED const ArraySeqpos *storeonline,
                           const ArraySeqpos *storeoffline)
{
  unsigned long ss;

  assert(storeoffline->nextfreeSeqpos == storeonline->nextfreeSeqpos);
  if (storeoffline->nextfreeSeqpos > 1UL)
  {
    qsort(storeoffline->spaceSeqpos,(size_t) storeoffline->nextfreeSeqpos,
          sizeof (Seqpos),
          cmpdescend);
  }
  for (ss=0; ss < storeoffline->nextfreeSeqpos; ss++)
  {
    assert(storeoffline->spaceSeqpos != NULL &&
           storeonline->spaceSeqpos != NULL);
    assert(storeoffline->spaceSeqpos[ss] == storeonline->spaceSeqpos[ss]);
  }
}

static void complementtag(Uchar *transformedtag,unsigned long taglen)
{
  unsigned long idx;

  for (idx = 0; idx < taglen; idx++)
  {
    transformedtag[taglen - 1 - idx] = COMPLEMENTBASE(transformedtag[idx]);
  }
}

/*
  XXX add forward and reverse match in one iteration
*/

int runtagerator(const TageratorOptions *tageratoroptions,Error *err)
{
  Suffixarray suffixarray;
  Seqpos totallength;
  SeqIterator *seqit = NULL;
  bool haserr = false;
  int retval, try;
  unsigned int demand = SARR_SUFTAB | SARR_ESQTAB;
  Limdfsresources *limdfsresources = NULL;
  Myersonlineresources *mor = NULL;
  ArraySeqpos storeonline, storeoffline;

  if (mapsuffixarray(&suffixarray,
                     &totallength,
                     demand,
                     tageratoroptions->indexname,
                     NULL,
                     err) != 0)
  {
    haserr = true;
  }
  if (suffixarray.readmode != Forwardmode)
  {
    error_set(err,"can only process index in forward mode");
    haserr = true;
  }
  INITARRAY(&storeonline,Seqpos);
  INITARRAY(&storeoffline,Seqpos);
  if (!haserr)
  {
    unsigned long taglen;
    uint64_t tagnumber;
    unsigned int mapsize;
    const Uchar *symbolmap, *currenttag;
    Uchar transformedtag[MAXTAGSIZE];
    char *desc = NULL;
    void (*processmatch)(void *,Seqpos,Seqpos);
    void *processmatchinfoonline, *processmatchinfooffline;

    symbolmap = getsymbolmapAlphabet(suffixarray.alpha);
    mapsize = getmapsizeAlphabet(suffixarray.alpha);
    if (tageratoroptions->docompare)
    {
      processmatch = storematch;
      processmatchinfoonline = &storeonline;
      processmatchinfooffline = &storeoffline;
    } else
    {
      processmatch = showmatch;
      processmatchinfoonline = NULL;
      processmatchinfooffline = NULL;
    }
    if (tageratoroptions->online || tageratoroptions->docompare)
    {
      mor = newMyersonlineresources(mapsize,suffixarray.encseq,
                                    processmatch,
                                    processmatchinfoonline);
    }
    if (tageratoroptions->maxdistance > 0)
    {
      limdfsresources = newLimdfsresources(suffixarray.encseq,
                                           suffixarray.readmode,
                                           mapsize,
                                           suffixarray.suftab,
                                           processmatch,
                                           processmatchinfooffline);
    }
    seqit = seqiterator_new(tageratoroptions->tagfiles, NULL, true);
    for (tagnumber = 0; !haserr; tagnumber++)
    {
      retval = seqiterator_next(seqit, &currenttag, &taglen, &desc, err);
      if (retval != 1)
      {
        if (retval < 0)
        {
          ma_free(desc);
        }
        break;
      }
      if (dotransformtag(transformedtag,
                         symbolmap,
                         currenttag,
                         taglen,
                         tagnumber,
                         tageratoroptions->replacewildcard,
                         err) != 0)
      {
        haserr = true;
        ma_free(desc);
        break;
      }
      printf("# patternlength=%lu\n",taglen);
      printf("# maxdistance=%lu\n",tageratoroptions->maxdistance);
      printf("# tag=");
      showsymbolstringgeneric(stdout,suffixarray.alpha,transformedtag,taglen);
      printf("\n");
      storeoffline.nextfreeSeqpos = 0;
      storeonline.nextfreeSeqpos = 0;
      assert(taglen > tageratoroptions->maxdistance);
      for (try=0 ; try < 2; try++)
      {
        if ((try == 0 && tageratoroptions->fwdmatch) ||
            (try == 1 && tageratoroptions->rcmatch))
        {
          if (try == 1 && tageratoroptions->rcmatch)
          {
            complementtag(transformedtag,taglen);
          }
          performthesearch(tageratoroptions,
                           mor,
                           limdfsresources,
                           suffixarray.encseq,
                           suffixarray.readmode,
                           suffixarray.suftab,
                           totallength,
                           transformedtag,
                           taglen,
                           processmatch,
                           processmatchinfooffline,
                           (try == 1 && tageratoroptions->rcmatch)
                             ? true : false);
          if (tageratoroptions->docompare)
          {
            compareresults(&storeonline,&storeoffline);
          }
        }
      }
      ma_free(desc);
    }
  }
  FREEARRAY(&storeonline,Seqpos);
  FREEARRAY(&storeoffline,Seqpos);
  if (limdfsresources != NULL)
  {
    freeLimdfsresources(&limdfsresources);
  }
  if (mor != NULL)
  {
    freeMyersonlineresources(&mor);
  }
  seqiterator_delete(seqit);
  freesuffixarray(&suffixarray);
  return haserr ? -1 : 0;
}
