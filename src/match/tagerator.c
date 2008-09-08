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
#include "core/unused_api.h"
#include "core/strarray.h"
#include "core/ma.h"
#include "core/error.h"
#include "core/fileutils.h"
#include "core/seqiterator.h"
#include "core/arraydef.h"
#include "tagerator.h"
#include "sarr-def.h"
#include "intbits.h"
#include "alphadef.h"
#include "myersapm.h"
#include "eis-voiditf.h"
#include "format64.h"
#include "idx-limdfs.h"
#include "mssufpat.h"
#include "apmoveridx.h"
#include "dist-short.h"
#include "stamp.h"

#include "echoseq.pr"
#include "esa-map.pr"

#define MAXTAGSIZE INTWORDSIZE

typedef struct
{
  Seqpos dbstartpos,
         matchlength;
  bool rcmatch;
} Simplematch;

typedef struct
{
  Uchar transformedtag[MAXTAGSIZE];
  unsigned long taglen;
  bool rcdir;
} Tagwithlength;

typedef struct
{
  const TageratorOptions *tageratoroptions;
  unsigned int alphasize;
  const Uchar *tagptr;
  const Alphabet *alpha;
  unsigned long *eqsvector;
} Showmatchinfo;

static void showmatch(void *processinfo,
                      bool rcmatch,
                      Seqpos dbstartpos,
                      Seqpos dblen,
                      const Uchar *dbsubstring,
                      unsigned long pprefixlen)
{
  Showmatchinfo *showmatchinfo = (Showmatchinfo *) processinfo;

  printf(FormatSeqpos,PRINTSeqposcast(dblen));
  printf(" %c " FormatSeqpos,rcmatch ? '-' : '+',PRINTSeqposcast(dbstartpos));
  if (showmatchinfo->tageratoroptions != NULL &&
      showmatchinfo->tageratoroptions->maxintervalwidth > 0)
  {
    printf(" ");
    printfsymbolstring(showmatchinfo->alpha,dbsubstring,(unsigned long) dblen);
    if (showmatchinfo->tageratoroptions->skpp)
    {
      unsigned long suffixlength
        = reversesuffixmatch(showmatchinfo->eqsvector,
                             showmatchinfo->alphasize,
                             dbsubstring,
                             (unsigned long) dblen,
                             showmatchinfo->tagptr,
                             pprefixlen,
                             (unsigned long) showmatchinfo->tageratoroptions
                                                         ->maxdistance);
      assert(pprefixlen >= suffixlength);
      printf(" %lu %lu ",suffixlength,pprefixlen - suffixlength);
      printfsymbolstring(NULL,showmatchinfo->tagptr +
                              (pprefixlen - suffixlength),
                              suffixlength);
    } else
    {
      printf(" %lu 0 ",pprefixlen);
      printfsymbolstring(NULL,showmatchinfo->tagptr, pprefixlen);
    }
  }
  printf("\n");
}

DECLAREARRAYSTRUCT(Simplematch);

static void storematch(void *processinfo,
                       bool rcmatch,
                       Seqpos dbstartpos,
                       Seqpos dblen,
                       GT_UNUSED const Uchar *dbsubstring,
                       GT_UNUSED unsigned long pprefixlen)
{
  ArraySimplematch *storetab = (ArraySimplematch *) processinfo;
  Simplematch *match;

  GETNEXTFREEINARRAY(match,storetab,Simplematch,32);
  match->dbstartpos = dbstartpos;
  match->matchlength = dblen;
  match->rcmatch = rcmatch;
}

static void checkmstats(void *processinfo,
                        const void *patterninfo,
                        unsigned long patternstartpos,
                        unsigned long mstatlength,
                        Seqpos leftbound,
                        Seqpos rightbound)
{
  unsigned long realmstatlength;
  Tagwithlength *twl = (Tagwithlength *) patterninfo;

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
  if (intervalwidthleq((const Limdfsresources *) processinfo,leftbound,
                       rightbound))
  {
    Uchar cc;
    Seqpos *sptr, witnessposition;
    unsigned long idx;
    ArraySeqpos *mstatspos = fromitv2sortedmatchpositions(
                                  (Limdfsresources *) processinfo,
                                  twl->rcdir,
                                  leftbound,
                                  rightbound,
                                  mstatlength);
    for (sptr = mstatspos->spaceSeqpos; sptr < mstatspos->spaceSeqpos +
                                               mstatspos->nextfreeSeqpos;
         sptr++)
    {
      witnessposition = *sptr;
      for (idx = patternstartpos; idx < patternstartpos + mstatlength; idx++)
      {
        cc = limdfsgetencodedchar((const Limdfsresources *) processinfo,
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
                          (unsigned long)
                          (witnessposition+idx-patternstartpos));
          exit(EXIT_FAILURE);
        }
      }
    }
  }
}

static void showmstats(void *processinfo,
                       const void *patterninfo,
                       GT_UNUSED unsigned long patternstartpos,
                       unsigned long mstatlength,
                       Seqpos leftbound,
                       Seqpos rightbound)
{
  Tagwithlength *twl = (Tagwithlength *) patterninfo;

  printf("%lu %c",mstatlength,twl->rcdir ? '-' : '+');
  if (intervalwidthleq((const Limdfsresources *) processinfo,leftbound,
                       rightbound))
  {
    unsigned long idx;
    ArraySeqpos *mstatspos = fromitv2sortedmatchpositions(
                                  (Limdfsresources *) processinfo,
                                  twl->rcdir,
                                  leftbound,
                                  rightbound,
                                  mstatlength);
    for (idx = 0; idx<mstatspos->nextfreeSeqpos; idx++)
    {
      printf(" " FormatSeqpos,PRINTSeqposcast(mstatspos->spaceSeqpos[idx]));
    }
  }
  printf("\n");
}

static int cmpdescend(const void *a,const void *b)
{
  Simplematch *valuea = (Simplematch *) a;
  Simplematch *valueb = (Simplematch *) b;

  if (!valuea->rcmatch && valueb->rcmatch)
  {
    return -1;
  }
  if (valuea->rcmatch && !valueb->rcmatch)
  {
    return 1;
  }
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
                          GT_Error *err)
{
  unsigned long idx;
  Uchar charcode;

  if (taglen > (unsigned long) MAXTAGSIZE)
  {
    gt_error_set(err,"tag \"%*.*s\" of length %lu; "
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
      gt_error_set(err,"undefined character '%c' in tag number " Formatuint64_t,
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
        gt_error_set(err,"wildcard in tag number " Formatuint64_t,
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
                                 bool rcmatch,
                                 Processmatch processmatch,
                                 void *processmatchinfooffline)
{
  if (tageratoroptions->online || (tageratoroptions->maxdistance >= 0 &&
                                   tageratoroptions->docompare))
  {
    edistmyersbitvectorAPM(mor,
                           transformedtag,
                           taglen,
                           rcmatch,
                           (unsigned long) tageratoroptions->maxdistance);
  }
  if (!tageratoroptions->online || tageratoroptions->docompare)
  {
    if (tageratoroptions->maxdistance == 0)
    {
      indexbasedexactpatternmatching(limdfsresources,
                                     rcmatch,
                                     transformedtag,
                                     taglen,
                                     processmatch,
                                     processmatchinfooffline);
    } else
    {
      indexbasedapproxpatternmatching(limdfsresources,
                                      rcmatch,
                                      transformedtag,
                                      taglen,
                                      (tageratoroptions->maxdistance < 0)
                                        ?  0
                                        : (unsigned long)
                                          tageratoroptions->maxdistance,
                                      tageratoroptions->maxintervalwidth,
                                      tageratoroptions->skpp,
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
    if (storeonline->spaceSimplematch[ss].rcmatch &&
        !storeoffline->spaceSimplematch[ss].rcmatch)
    {
      fprintf(stderr,"rcmatch: storeonline[%lu] = p != d "
                     "= storeoffline[%lu]\n",ss,ss);
      exit(EXIT_FAILURE);
    }
    if (!storeonline->spaceSimplematch[ss].rcmatch &&
        storeoffline->spaceSimplematch[ss].rcmatch)
    {
      fprintf(stderr,"rcmatch: storeonline[%lu] = d != p "
                     "= storeoffline[%lu]\n",ss,ss);
      exit(EXIT_FAILURE);
    }
    if (storeonline->spaceSimplematch[ss].matchlength !=
        storeoffline->spaceSimplematch[ss].matchlength)
    {
      fprintf(stderr,"matchlength: storeonline[%lu] = " FormatSeqpos
                     " != " FormatSeqpos "= storeoffline[%lu]\n",
                     ss,
                     PRINTSeqposcast(storeonline->spaceSimplematch[ss].
                                     matchlength),
                     PRINTSeqposcast(storeoffline->spaceSimplematch[ss].
                                     matchlength),
                     ss);
      exit(EXIT_FAILURE);
    }
    if (storeonline->spaceSimplematch[ss].dbstartpos !=
        storeoffline->spaceSimplematch[ss].dbstartpos)
    {
      fprintf(stderr,"dbstartpos: storeonline[%lu] = " FormatSeqpos
                     " != " FormatSeqpos "= storeoffline[%lu]\n",
                     ss,
                     PRINTSeqposcast(storeonline->spaceSimplematch[ss].
                                     dbstartpos),
                     PRINTSeqposcast(storeoffline->spaceSimplematch[ss].
                                     dbstartpos),
                     ss);
      exit(EXIT_FAILURE);
    }
  }
}

static void reversecomplementtag(Uchar *transformedtag,unsigned long taglen)
{
  Uchar tmp, *frontptr, *backptr;

  for (frontptr = transformedtag, backptr = transformedtag + taglen - 1;
       frontptr < backptr; frontptr++, backptr--)
  {
    tmp = *frontptr;
    *frontptr = COMPLEMENTBASE(*backptr);
    *backptr = COMPLEMENTBASE(tmp);
  }
}

int runtagerator(const TageratorOptions *tageratoroptions,GT_Error *err)
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
  Showmatchinfo showmatchinfo;

  if (tageratoroptions->maxdistance >= 0)
  {
    dfst = apm_AbstractDfstransformer();
  } else
  {
    dfst = pms_AbstractDfstransformer();
  }
  if (gt_str_length(tageratoroptions->esaindexname) > 0)
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
      gt_error_set(err,"using option -esa you can only process index "
                    "in forward mode");
      haserr = true;
    } else
    {
      if (!withesa && suffixarray.readmode != Reversemode)
      {
        gt_error_set(err,"with option -pck you can only process index "
                      "in reverse mode");
        haserr = true;
      }
    }
  }
  if (!haserr && gt_str_length(tageratoroptions->pckindexname) > 0)
  {
    packedindex = loadvoidBWTSeqForSA(tageratoroptions->pckindexname,
                                      &suffixarray,
                                      totallength, true, err);
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
    const Matchbound **mbtab;
    unsigned int maxdepth;
    unsigned long maxpathlength;
    Processmatch processmatch;
    void *processmatchinfoonline, *processmatchinfooffline;

    symbolmap = getsymbolmapAlphabet(suffixarray.alpha);
    mapsize = getmapsizeAlphabet(suffixarray.alpha);
    if (tageratoroptions->docompare)
    {
      processmatch = storematch;
      processmatchinfoonline = &storeonline;
      processmatchinfooffline = &storeoffline;
      showmatchinfo.eqsvector = NULL;
    } else
    {
      processmatch = showmatch;
      processmatchinfoonline = NULL;
      showmatchinfo.tageratoroptions = tageratoroptions;
      showmatchinfo.alphasize = (unsigned int) (mapsize-1);
      showmatchinfo.tagptr = &twl.transformedtag[0];
      showmatchinfo.alpha = suffixarray.alpha;
      showmatchinfo.eqsvector = gt_malloc(sizeof(*showmatchinfo.eqsvector) *
                                          showmatchinfo.alphasize);
      processmatchinfooffline = &showmatchinfo;
    }
    if (tageratoroptions->online || tageratoroptions->docompare)
    {
      assert(suffixarray.encseq != NULL);
      mor = newMyersonlineresources(mapsize,
                                    tageratoroptions->nowildcards,
                                    suffixarray.encseq,
                                    processmatch,
                                    processmatchinfoonline);
    }
    if (withesa)
    {
      mbtab = NULL;
      maxdepth = 0;
    } else
    {
      mbtab = bwtseq2mbtab(packedindex);
      maxdepth = bwtseq2maxdepth(packedindex);
      if (tageratoroptions->userdefinedmaxdepth >= 0 &&
          maxdepth > (unsigned int) tageratoroptions->userdefinedmaxdepth)
      {
        maxdepth = (unsigned int) tageratoroptions->userdefinedmaxdepth;
      }
    }
    if (tageratoroptions->maxdistance >= 0)
    {
      maxpathlength = (unsigned long) (1+ MAXTAGSIZE +
                                       tageratoroptions->maxdistance);
    } else
    {
      maxpathlength = (unsigned long) (1+MAXTAGSIZE);
    }
    limdfsresources = newLimdfsresources(withesa ? &suffixarray : packedindex,
                                         mbtab,
                                         maxdepth,
                                         suffixarray.encseq,
                                         withesa,
                                         tageratoroptions->nowildcards,
                                         tageratoroptions->maxintervalwidth,
                                         mapsize,
                                         totallength,
                                         maxpathlength,
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
          gt_free(desc);
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
        gt_free(desc);
        break;
      }
      twl.rcdir = false;
      printf("# %lu ",twl.taglen);
      fprintfsymbolstring(stdout,suffixarray.alpha,twl.transformedtag,
                          twl.taglen);
      printf("\n");
      storeoffline.nextfreeSimplematch = 0;
      storeonline.nextfreeSimplematch = 0;
      if (tageratoroptions->maxdistance > 0 &&
          twl.taglen <= (unsigned long) tageratoroptions->maxdistance)
      {
        gt_error_set(err,"tag \"%*.*s\" of length %lu; "
                  "tags must be longer than the allowed number of errors "
                  "(which is %ld)",
                   (int) twl.taglen,(int) twl.taglen,currenttag,twl.taglen,
                   tageratoroptions->maxdistance);
        haserr = true;
        gt_free(desc);
        break;
      }
      assert(tageratoroptions->maxdistance < 0 ||
             twl.taglen > (unsigned long) tageratoroptions->maxdistance);
      for (try=0 ; try < 2; try++)
      {
        if ((try == 0 && !tageratoroptions->nofwdmatch) ||
            (try == 1 && !tageratoroptions->norcmatch))
        {
          if (try == 1 && !tageratoroptions->norcmatch)
          {
            reversecomplementtag(twl.transformedtag,twl.taglen);
            twl.rcdir = true;
          }
          performpatternsearch(dfst,
                               tageratoroptions,
                               mor,
                               limdfsresources,
                               twl.transformedtag,
                               twl.taglen,
                               twl.rcdir,
                               processmatch,
                               processmatchinfooffline);
          if (tageratoroptions->docompare)
          {
            compareresults(&storeonline,&storeoffline);
          }
        }
      }
      gt_free(desc);
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
  gt_free(showmatchinfo.eqsvector);
  if (packedindex != NULL)
  {
    deletevoidBWTSeq(packedindex);
  }
  return haserr ? -1 : 0;
}
