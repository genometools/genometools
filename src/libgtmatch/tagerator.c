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
#include "intbits.h"
#include "alphadef.h"
#include "esa-myersapm.h"
#include "eis-voiditf.h"
#include "format64.h"
#include "idx-limdfs.h"
#include "apmoveridx.h"
#include "mssufpat.h"
#include "stamp.h"

#include "echoseq.pr"
#include "esa-map.pr"

#define MAXTAGSIZE INTWORDSIZE

typedef struct
{
  Seqpos dbstartpos, matchlength;
} Simplematch;

typedef struct
{
  Uchar transformedtag[MAXTAGSIZE];
  unsigned long taglen;
} Tagwithlength;

static Seqpos convertstartpos(bool withesa,Seqpos totallength,Seqpos startpos)
{
  if (withesa)
  {
    return startpos;
  }
  assert(totallength >= startpos);
  return totallength - startpos;
}

static void showmatch(void *processinfo,bool withesa,
                      Seqpos totallength,Seqpos startpos,Seqpos len,
                      unsigned long pprefixlen)
{
  printf(FormatSeqpos,PRINTSeqposcast(len));
  printf(" " FormatSeqpos,
         PRINTSeqposcast(convertstartpos(withesa,totallength,startpos)));
  if (processinfo != NULL && *((unsigned long *) processinfo) > 0)
  {
    printf(" %lu",pprefixlen);
  }
  printf("\n");
}

static void checkmstats(void *processinfo,
                        const void *patterninfo,
                        unsigned long patternstartpos,
                        unsigned long mstatlength,
                        Seqpos leftbound)
{
  unsigned long realmstatlength;
  Tagwithlength *twl = (Tagwithlength *) patterninfo;
  Seqpos witnessposition;
  unsigned long idx;

  realmstatlength = genericmstats((const Limdfsresources *) processinfo,
                                  &twl->transformedtag[patternstartpos],
                                  &twl->transformedtag[twl->taglen]);
  if (mstatlength != realmstatlength)
  {
    fprintf(stderr,"patternstartpos = %lu: mstatlength = %lu != %lu "
                   " = realmstatlength\n",
                    patternstartpos,mstatlength,realmstatlength);
    exit(EXIT_FAILURE);
  }
  witnessposition = bound2startpos((const Limdfsresources *) processinfo,
                                   leftbound,mstatlength);
  for (idx = patternstartpos; idx < patternstartpos + mstatlength; idx++)
  {
    Uchar cc = limdfsgetencodedchar((const Limdfsresources *) processinfo,
                                    witnessposition + idx - patternstartpos,
                                    Forwardmode);
    if (twl->transformedtag[idx] != cc)
    {
      fprintf(stderr,"patternstartpos = %lu: pattern[%lu] = %u != %u = "
                     "sequence[%lu]\n",
                      patternstartpos,
                      idx,
                      (unsigned int) twl->transformedtag[idx],
                      (unsigned int) cc,
                      (unsigned long) (witnessposition+idx-patternstartpos));
      exit(EXIT_FAILURE);
    }
  }
}

static void showmstats(void *processinfo,
                       UNUSED const void *patterninfo,
                       UNUSED unsigned long idx,
                       unsigned long value,
                       Seqpos leftbound)
{
  Seqpos witnessposition = bound2startpos((const Limdfsresources *) processinfo,
                                          leftbound,value);
  printf("%lu " FormatSeqpos "\n",value,PRINTSeqposcast(witnessposition));
}

DECLAREARRAYSTRUCT(Simplematch);

static void storematch(void *processinfo,
                       bool withesa,
                       Seqpos totallength,
                       Seqpos startpos,
                       Seqpos len,
                       UNUSED unsigned long pprefixlen)
{
  ArraySimplematch *storetab = (ArraySimplematch *) processinfo;
  Simplematch *match;

  GETNEXTFREEINARRAY(match,storetab,Simplematch,32);
  match->dbstartpos = convertstartpos(withesa,totallength,startpos);
  match->matchlength = len;
}

static int cmpdescend(const void *a,const void *b)
{
  Simplematch *valuea = (Simplematch *) a;
  Simplematch *valueb = (Simplematch *) b;

  if (valuea->dbstartpos < valueb->dbstartpos)
  {
    return 1;
  }
  if (valuea->dbstartpos > valueb->dbstartpos)
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
    error_set(err,"tag \"%*.*s\" of length %lu; "
                  "tags must not be longer than %lu",
                   (int) taglen,(int) taglen,currenttag,taglen,
                   (unsigned long) MAXTAGSIZE);
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

static void performpatternsearch(const AbstractDfstransformer *dfst,
                                 const TageratorOptions *tageratoroptions,
                                 Myersonlineresources *mor,
                                 Limdfsresources *limdfsresources,
                                 const Uchar *transformedtag,
                                 unsigned long taglen,
                                 void (*processmatch)(void *,bool,Seqpos,
                                                      Seqpos,Seqpos,
                                                      unsigned long),
                                 void *processmatchinfooffline,
                                 UNUSED bool rcmatch)
{
  if (tageratoroptions->online || (tageratoroptions->maxdistance >= 0 &&
                                   tageratoroptions->docompare))
  {
    edistmyersbitvectorAPM(mor,
                           transformedtag,
                           taglen,
                           (unsigned long) tageratoroptions->maxdistance);
  }
  if (!tageratoroptions->online || tageratoroptions->docompare)
  {
    if (tageratoroptions->maxdistance == 0)
    {
      indexbasedexactpatternmatching(limdfsresources,
                                     transformedtag,
                                     taglen,
                                     processmatch,
                                     processmatchinfooffline);
    } else
    {
      indexbasedapproxpatternmatching(limdfsresources,
                                      transformedtag,
                                      taglen,
                                      (tageratoroptions->maxdistance < 0)
                                        ?  0
                                        : (unsigned long)
                                          tageratoroptions->maxdistance,
                                      (Seqpos)
                                          tageratoroptions->maxintervalwidth,
                                      dfst);
    }
  }
}

static void compareresults(const ArraySimplematch *storeonline,
                           const ArraySimplematch *storeoffline)
{
  unsigned long ss;

  if (storeonline->nextfreeSimplematch != storeoffline->nextfreeSimplematch)
  {
    fprintf(stderr,"nextfreeSimplematch: storeonline = %lu != %lu "
                   "storeoffline\n",
                   storeonline->nextfreeSimplematch,
                   storeoffline->nextfreeSimplematch);
    exit(EXIT_FAILURE);
  }
  assert(storeonline->nextfreeSimplematch == storeoffline->nextfreeSimplematch);
  if (storeoffline->nextfreeSimplematch > 1UL)
  {
    qsort(storeoffline->spaceSimplematch,(size_t)
          storeoffline->nextfreeSimplematch,
          sizeof (Simplematch),
          cmpdescend);
  }
  for (ss=0; ss < storeoffline->nextfreeSimplematch; ss++)
  {
    assert(storeonline->spaceSimplematch != NULL &&
           storeoffline->spaceSimplematch != NULL);
    if (storeonline->spaceSimplematch[ss].matchlength !=
        storeoffline->spaceSimplematch[ss].matchlength)
    {
      fprintf(stderr,"matchlength: storeonline[%lu] = %lu != %lu "
                     "= storeoffline[%lu]\n",
                     ss,
                     (unsigned long)
                     storeonline->spaceSimplematch[ss].matchlength,
                     (unsigned long)
                     storeoffline->spaceSimplematch[ss].matchlength,
                     ss);
    }
    assert(storeoffline->spaceSimplematch[ss].dbstartpos ==
           storeonline->spaceSimplematch[ss].dbstartpos);
    assert(storeoffline->spaceSimplematch[ss].matchlength ==
           storeonline->spaceSimplematch[ss].matchlength);
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
  unsigned int demand;
  Limdfsresources *limdfsresources = NULL;
  Myersonlineresources *mor = NULL;
  ArraySimplematch storeonline, storeoffline;
  void *packedindex = NULL;
  bool withesa;
  const AbstractDfstransformer *dfst;

  if (tageratoroptions->maxdistance >= 0)
  {
    dfst = apm_AbstractDfstransformer();
  } else
  {
    dfst = pms_AbstractDfstransformer();
  }
  if (str_length(tageratoroptions->esaindexname) > 0)
  {
    demand = SARR_ESQTAB;
    if (!tageratoroptions->online)
    {
      demand |= SARR_SUFTAB;
    }
    withesa = true;
  } else
  {
    if (tageratoroptions->docompare || tageratoroptions->online)
    {
      demand = SARR_ESQTAB;
    } else
    {
      demand = 0;
    }
    withesa = false;
  }
  if (mapsuffixarray(&suffixarray,
                     &totallength,
                     demand,
                     withesa ? tageratoroptions->esaindexname
                             : tageratoroptions->pckindexname,
                     NULL,
                     err) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (withesa && suffixarray.readmode != Forwardmode)
    {
      error_set(err,"using option -esa you can only process index "
                    "in forward mode");
      haserr = true;
    } else
    {
      if (!withesa && suffixarray.readmode != Reversemode)
      {
        error_set(err,"with option -pck you can only process index "
                      "in reverse mode");
        haserr = true;
      }
    }
  }
  if (!haserr && str_length(tageratoroptions->pckindexname) > 0)
  {
    packedindex = loadvoidBWTSeqForSA(tageratoroptions->pckindexname,
                                      &suffixarray,
                                      totallength, err);
    if (packedindex == NULL)
    {
      haserr = true;
    }
  }
  INITARRAY(&storeonline,Simplematch);
  INITARRAY(&storeoffline,Simplematch);
  if (!haserr)
  {
    Tagwithlength twl;
    uint64_t tagnumber;
    unsigned int mapsize;
    const Uchar *symbolmap, *currenttag;
    char *desc = NULL;
    void (*processmatch)(void *,bool,Seqpos,Seqpos,Seqpos,unsigned long);
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
      processmatchinfooffline = (void *) &tageratoroptions->maxintervalwidth;
    }
    if (tageratoroptions->online || tageratoroptions->docompare)
    {
      assert(suffixarray.encseq != NULL);
      mor = newMyersonlineresources(mapsize,
                                    tageratoroptions->nospecials,
                                    suffixarray.encseq,
                                    processmatch,
                                    processmatchinfoonline);
    }
    limdfsresources = newLimdfsresources(withesa ? &suffixarray : packedindex,
                                         suffixarray.encseq,
                                         withesa,
                                         tageratoroptions->nospecials,
                                         mapsize,
                                         totallength,
                                         processmatch,
                                         processmatchinfooffline,
                                         tageratoroptions->docompare
                                           ? checkmstats
                                           : showmstats,
                                         &twl, /* refer to uninit structure */
                                         dfst);
    seqit = seqiterator_new(tageratoroptions->tagfiles, NULL, true);
    for (tagnumber = 0; !haserr; tagnumber++)
    {
      retval = seqiterator_next(seqit, &currenttag, &twl.taglen, &desc, err);
      if (retval != 1)
      {
        if (retval < 0)
        {
          ma_free(desc);
        }
        break;
      }
      if (dotransformtag(&twl.transformedtag[0],
                         symbolmap,
                         currenttag,
                         twl.taglen,
                         tagnumber,
                         tageratoroptions->replacewildcard,
                         err) != 0)
      {
        haserr = true;
        ma_free(desc);
        break;
      }
      printf("# %lu ",twl.taglen);
      fprintfsymbolstring(stdout,suffixarray.alpha,twl.transformedtag,
                          twl.taglen);
      printf("\n");
      storeoffline.nextfreeSimplematch = 0;
      storeonline.nextfreeSimplematch = 0;
      assert(tageratoroptions->maxdistance < 0 ||
             twl.taglen > (unsigned long) tageratoroptions->maxdistance);
      for (try=0 ; try < 2; try++)
      {
        if ((try == 0 && tageratoroptions->fwdmatch) ||
            (try == 1 && tageratoroptions->rcmatch))
        {
          if (try == 1 && tageratoroptions->rcmatch)
          {
            complementtag(twl.transformedtag,twl.taglen);
          }
          performpatternsearch(dfst,
                               tageratoroptions,
                               mor,
                               limdfsresources,
                               twl.transformedtag,
                               twl.taglen,
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
  FREEARRAY(&storeonline,Simplematch);
  FREEARRAY(&storeoffline,Simplematch);
  if (limdfsresources != NULL)
  {
    freeLimdfsresources(&limdfsresources,dfst);
  }
  if (mor != NULL)
  {
    freeMyersonlineresources(&mor);
  }
  seqiterator_delete(seqit);
  freesuffixarray(&suffixarray);
  if (packedindex != NULL)
  {
    deletevoidBWTSeq(packedindex);
  }
  return haserr ? -1 : 0;
}
