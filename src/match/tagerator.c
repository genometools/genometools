/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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
#include "core/str_array.h"
#include "core/ma.h"
#include "core/error.h"
#include "core/fileutils.h"
#include "core/seqiterator.h"
#include "core/arraydef.h"
#include "revcompl.h"
#include "sarr-def.h"
#include "intbits.h"
#include "alphadef.h"
#include "myersapm.h"
#include "format64.h"
#include "idx-limdfs.h"
#include "mssufpat.h"
#include "apmeoveridx.h"
#include "dist-short.h"
#include "stamp.h"
#include "tagerator.h"
#include "esa-map.h"
#include "echoseq.h"

#define MAXTAGSIZE INTWORDSIZE

#define ISRCDIR(TWL)  (((TWL)->tagptr == (TWL)->transformedtag)\
                        ? false\
                        : true)

typedef struct
{
  Seqpos dbstartpos,
         matchlength;
  bool rcmatch;
} Simplematch;

typedef struct
{
  const Uchar *tagptr;
  Uchar transformedtag[MAXTAGSIZE],
        rctransformedtag[MAXTAGSIZE];
  unsigned long taglen;
} Tagwithlength;

typedef struct
{
  const TageratorOptions *tageratoroptions;
  unsigned int alphasize;
  const Uchar *tagptr;
  const SfxAlphabet *alpha;
  unsigned long *eqsvector;
  const Tagwithlength *twlptr;
} Showmatchinfo;

#define ADDTABULATOR\
        if (firstitem)\
        {\
          firstitem = false;\
        } else\
        {\
          (void) putchar('\t');\
        }

static void showmatch(void *processinfo,
                      const GtMatch *match)
{
  Showmatchinfo *showmatchinfo = (Showmatchinfo *) processinfo;
  bool firstitem = true;

  gt_assert(showmatchinfo->tageratoroptions != NULL);
  if (showmatchinfo->tageratoroptions->outputmode & TAGOUT_DBLENGTH)
  {
    printf(FormatSeqpos,PRINTSeqposcast(match->dblen));
    firstitem = false;
  }
  if (showmatchinfo->tageratoroptions->outputmode & TAGOUT_DBSTARTPOS)
  {
    ADDTABULATOR;
    printf(FormatSeqpos,PRINTSeqposcast(match->dbstartpos));
  }
  if (showmatchinfo->tageratoroptions->outputmode & TAGOUT_DBSEQUENCE)
  {
    ADDTABULATOR;
    gt_assert(match->dbsubstring != NULL);
    printfsymbolstring(showmatchinfo->alpha,match->dbsubstring,
                       (unsigned long) match->dblen);
  }
  if (showmatchinfo->tageratoroptions->outputmode & TAGOUT_STRAND)
  {
    ADDTABULATOR;
    printf("%c",ISRCDIR(showmatchinfo->twlptr) ? '-' : '+');
  }
  if (showmatchinfo->tageratoroptions->outputmode & TAGOUT_EDIST)
  {
    ADDTABULATOR;
    printf("%lu",match->distance);
  }
  if (showmatchinfo->tageratoroptions->maxintervalwidth > 0)
  {
    if (showmatchinfo->tageratoroptions->skpp)
    {
      if (showmatchinfo->tageratoroptions->outputmode &
          (TAGOUT_TAGSTARTPOS | TAGOUT_TAGLENGTH | TAGOUT_TAGSUFFIXSEQ))
      {
        unsigned long suffixlength
          = reversesuffixmatch(showmatchinfo->eqsvector,
                               showmatchinfo->alphasize,
                               match->dbsubstring,
                               (unsigned long) match->dblen,
                               showmatchinfo->tagptr,
                               match->pprefixlen,
                               (unsigned long) showmatchinfo->tageratoroptions->
                                                        userdefinedmaxdistance);
        gt_assert(match->pprefixlen >= suffixlength);
        if (showmatchinfo->tageratoroptions->outputmode & TAGOUT_TAGSTARTPOS)
        {
          ADDTABULATOR;
          printf("%lu",match->pprefixlen - suffixlength);
        }
        if (showmatchinfo->tageratoroptions->outputmode & TAGOUT_TAGLENGTH)
        {
          ADDTABULATOR;
          printf("%lu",suffixlength);
        }
        if (showmatchinfo->tageratoroptions->outputmode & TAGOUT_TAGSUFFIXSEQ)
        {
          ADDTABULATOR;
          printfsymbolstring(NULL,showmatchinfo->tagptr +
                                  (match->pprefixlen - suffixlength),
                                  suffixlength);
        }
      }
    } else
    {
      if (showmatchinfo->tageratoroptions->outputmode & TAGOUT_TAGSTARTPOS)
      {
        ADDTABULATOR;
        printf("0");
      }
      if (showmatchinfo->tageratoroptions->outputmode & TAGOUT_TAGLENGTH)
      {
        ADDTABULATOR;
        printf("%lu",match->pprefixlen);
      }
      if (showmatchinfo->tageratoroptions->outputmode & TAGOUT_TAGSUFFIXSEQ)
      {
        ADDTABULATOR;
        printfsymbolstring(NULL,showmatchinfo->tagptr, match->pprefixlen);
      }
    }
  }
  if (!firstitem)
  {
    printf("\n");
  }
}

typedef struct
{
  Simplematch *spaceSimplematch;
  unsigned long nextfreeSimplematch, allocatedSimplematch;
  const Tagwithlength *twlptr;
} ArraySimplematch;

static void storematch(void *processinfo,const GtMatch *match)
{
  ArraySimplematch *storetab = (ArraySimplematch *) processinfo;
  Simplematch *simplematch;

  GETNEXTFREEINARRAY(simplematch,storetab,Simplematch,32);
  simplematch->dbstartpos = match->dbstartpos;
  simplematch->matchlength = match->dblen;
  simplematch->rcmatch = ISRCDIR(storetab->twlptr);
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
                                  twl->tagptr + patternstartpos,
                                  twl->tagptr + twl->taglen);
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
        if (twl->tagptr[idx] != cc)
        {
          fprintf(stderr,"patternstartpos = %lu: pattern[%lu] = %u != %u = "
                         "sequence[%lu]\n",
                          patternstartpos,
                          idx,
                          (unsigned int) twl->tagptr[idx],
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

  printf("%lu %c",mstatlength,ISRCDIR(twl) ? '-' : '+');
  if (intervalwidthleq((const Limdfsresources *) processinfo,leftbound,
                       rightbound))
  {
    unsigned long idx;
    ArraySeqpos *mstatspos = fromitv2sortedmatchpositions(
                                  (Limdfsresources *) processinfo,
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
                          GtError *err)
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
        charcode = 0;
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

static bool performpatternsearch(const AbstractDfstransformer *dfst,
                                 bool domstats,
                                 unsigned long maxdistance,
                                 bool online,
                                 bool docompare,
                                 unsigned long maxintervalwidth,
                                 bool skpp,
                                 Myersonlineresources *mor,
                                 Limdfsresources *limdfsresources,
                                 const Uchar *tagptr,
                                 unsigned long taglen)
{
  if (online || (!domstats && docompare))
  {
    gt_assert(mor != NULL);
    edistmyersbitvectorAPM(mor,tagptr,taglen,maxdistance);
  }
  if (!online || docompare)
  {
    if (domstats)
    {
      indexbasedmstats(limdfsresources,tagptr,taglen,dfst);
      return false;
    }
    if (maxdistance == 0)
    {
      return indexbasedexactpatternmatching(limdfsresources,tagptr,taglen);
    } else
    {
      return indexbasedapproxpatternmatching(limdfsresources,tagptr,taglen,
                                             maxdistance,
                                             maxintervalwidth,
                                             skpp,
                                             dfst);
    }
  }
  return false;
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
  gt_assert(storeonline->nextfreeSimplematch ==
            storeoffline->nextfreeSimplematch);
  if (storeoffline->nextfreeSimplematch > 1UL)
  {
    qsort(storeoffline->spaceSimplematch,(size_t)
          storeoffline->nextfreeSimplematch,
          sizeof (Simplematch),
          cmpdescend);
  }
  for (ss=0; ss < storeoffline->nextfreeSimplematch; ss++)
  {
    gt_assert(storeonline->spaceSimplematch != NULL &&
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

static void searchoverstrands(const TageratorOptions *tageratoroptions,
                              Tagwithlength *twl,
                              const AbstractDfstransformer *dfst,
                              Myersonlineresources *mor,
                              Limdfsresources *limdfsresources,
                              Showmatchinfo *showmatchinfo,
                              ArraySimplematch *storeonline,
                              ArraySimplematch *storeoffline)
{
  int try;
  bool domstats, matchfound;
  unsigned long maxdistance, mindistance, distance;

  if (tageratoroptions->userdefinedmaxdistance < 0)
  {
    domstats = true;
    mindistance = maxdistance = 0;
  } else
  {
    domstats = false;
    gt_assert(tageratoroptions->userdefinedmaxdistance >= 0);
    maxdistance = (unsigned long) tageratoroptions->userdefinedmaxdistance;
    if (tageratoroptions->best)
    {
      mindistance = 0;
    } else
    {
      mindistance = maxdistance;
    }
  }
  matchfound = false;
  for (distance = mindistance; distance <= maxdistance; distance++)
  {
    showmatchinfo->tagptr = twl->tagptr = twl->transformedtag;
    for (try=0 ; try < 2; try++)
    {
      if ((try == 0 && !tageratoroptions->nofwdmatch) ||
          (try == 1 && !tageratoroptions->norcmatch))
      {
        if (try == 1 && !tageratoroptions->norcmatch)
        {
          showmatchinfo->tagptr = twl->tagptr = twl->rctransformedtag;
        }
        if (performpatternsearch(dfst,
                                 domstats,
                                 distance,
                                 tageratoroptions->online,
                                 tageratoroptions->docompare,
                                 tageratoroptions->maxintervalwidth,
                                 tageratoroptions->skpp,
                                 mor,
                                 limdfsresources,
                                 twl->tagptr,
                                 twl->taglen) && !matchfound)
        {
          matchfound = true;
        }
        if (tageratoroptions->docompare)
        {
          compareresults(storeonline,storeoffline);
        }
      }
    }
    if (tageratoroptions->best && matchfound)
    {
      break;
    }
  }
}

int runtagerator(const TageratorOptions *tageratoroptions,GtError *err)
{
  GtSeqIterator *seqit = NULL;
  bool haserr = false;
  int retval;
  Myersonlineresources *mor = NULL;
  ArraySimplematch storeonline, storeoffline;
  const AbstractDfstransformer *dfst;
  Genericindex *genericindex = NULL;
  const Encodedsequence *encseq = NULL;
  Verboseinfo *verboseinfo;

  verboseinfo = newverboseinfo(tageratoroptions->verbose);
  if (tageratoroptions->userdefinedmaxdistance >= 0)
  {
    dfst = apme_AbstractDfstransformer();
  } else
  {
    dfst = pms_AbstractDfstransformer();
  }
  if (tageratoroptions->online)
  {
    encseq = mapencodedsequence (true,
                                 tageratoroptions->indexname,
                                 true,
                                 false,
                                 false,
                                 verboseinfo,
                                 err);
    if (encseq == NULL)
    {
      haserr = true;
    }
  } else
  {
    genericindex = genericindex_new(tageratoroptions->indexname,
                                    tageratoroptions->withesa,
                                    tageratoroptions->docompare,
                                    tageratoroptions->userdefinedmaxdepth,
                                    verboseinfo,
                                    err);
    if (genericindex == NULL)
    {
      haserr = true;
    } else
    {
      encseq = genericindex_getencseq(genericindex);
    }
  }
  INITARRAY(&storeonline,Simplematch);
  INITARRAY(&storeoffline,Simplematch);
  if (!haserr)
  {
    Tagwithlength twl;
    uint64_t tagnumber;
    unsigned int numofchars;
    const Uchar *symbolmap, *currenttag;
    char *desc = NULL;
    Processmatch processmatch;
    Showmatchinfo showmatchinfo;
    void *processmatchinfoonline, *processmatchinfooffline;
    Limdfsresources *limdfsresources = NULL;
    const SfxAlphabet *alpha;

    storeonline.twlptr = storeoffline.twlptr = &twl;
    alpha = getencseqAlphabet(encseq);
    symbolmap = getsymbolmapAlphabet(alpha);
    numofchars = getnumofcharsAlphabet(alpha);
    if (tageratoroptions->docompare)
    {
      processmatch = storematch;
      processmatchinfoonline = &storeonline;
      processmatchinfooffline = &storeoffline;
      showmatchinfo.eqsvector = NULL;
    } else
    {
      processmatch = showmatch;
      showmatchinfo.twlptr = &twl;
      showmatchinfo.tageratoroptions = tageratoroptions;
      showmatchinfo.alphasize = (unsigned int) numofchars;
      showmatchinfo.alpha = alpha;
      showmatchinfo.eqsvector = gt_malloc(sizeof(*showmatchinfo.eqsvector) *
                                          showmatchinfo.alphasize);
      processmatchinfooffline = &showmatchinfo;
      processmatchinfoonline = &showmatchinfo;
    }
    if (tageratoroptions->online || tageratoroptions->docompare)
    {
      gt_assert(encseq != NULL);
      mor = newMyersonlineresources(numofchars,
                                    tageratoroptions->nowildcards,
                                    encseq,
                                    processmatch,
                                    processmatchinfoonline);
    }
    if (!tageratoroptions->online || tageratoroptions->docompare)
    {
      unsigned long maxpathlength;

      if (tageratoroptions->userdefinedmaxdistance >= 0)
      {
        maxpathlength = (unsigned long) (1+ MAXTAGSIZE +
                                         tageratoroptions->
                                         userdefinedmaxdistance);
      } else
      {
        maxpathlength = (unsigned long) (1+MAXTAGSIZE);
      }
      limdfsresources = newLimdfsresources(genericindex,
                                           tageratoroptions->nowildcards,
                                           tageratoroptions->maxintervalwidth,
                                           maxpathlength,
                                           false, /* keepexpandedonstack */
                                           processmatch,
                                           processmatchinfooffline,
                                           tageratoroptions->docompare
                                             ? checkmstats
                                             : showmstats,
                                           &twl, /* refer to uninit structure */
                                           dfst);
    }
    printf("# for each match show: ");
    getsetargmodekeywords(tageratoroptions->modedesc,
                          tageratoroptions->numberofmodedescentries,
                          tageratoroptions->outputmode);
    seqit = gt_seqiterator_new(tageratoroptions->tagfiles, NULL, true);
    for (tagnumber = 0; !haserr; tagnumber++)
    {
      retval = gt_seqiterator_next(seqit, &currenttag, &twl.taglen, &desc, err);
      if (retval != 1)
      {
        if (retval < 0)
        {
          gt_free(desc);
        }
        break;
      }
      if (dotransformtag(twl.transformedtag,
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
      copy_reversecomplement(twl.rctransformedtag,twl.transformedtag,
                             twl.taglen);
      twl.tagptr = twl.transformedtag;
      printf("#");
      if (tageratoroptions->outputmode & TAGOUT_TAGNUM)
      {
        printf("\t" Formatuint64_t,PRINTuint64_tcast(tagnumber));
      }
      if (tageratoroptions->outputmode & TAGOUT_TAGSEQ)
      {
        printf("\t%lu\t",twl.taglen);
        fprintfsymbolstring(stdout,alpha,twl.transformedtag,twl.taglen);
      }
      printf("\n");
      storeoffline.nextfreeSimplematch = 0;
      storeonline.nextfreeSimplematch = 0;
      if (tageratoroptions->userdefinedmaxdistance > 0 &&
          twl.taglen <= (unsigned long)
                        tageratoroptions->userdefinedmaxdistance)
      {
        gt_error_set(err,"tag \"%*.*s\" of length %lu; "
                     "tags must be longer than the allowed number of errors "
                     "(which is %ld)",
                     (int) twl.taglen,
                     (int) twl.taglen,currenttag,
                     twl.taglen,
                     tageratoroptions->userdefinedmaxdistance);
        haserr = true;
        gt_free(desc);
        break;
      }
      gt_assert(tageratoroptions->userdefinedmaxdistance < 0 ||
                twl.taglen > (unsigned long)
                             tageratoroptions->userdefinedmaxdistance);
      searchoverstrands(tageratoroptions,
                        &twl,
                        dfst,
                        mor,
                        limdfsresources,
                        &showmatchinfo,
                        &storeonline,
                        &storeoffline);
      gt_free(desc);
    }
    gt_free(showmatchinfo.eqsvector);
    if (limdfsresources != NULL)
    {
      freeLimdfsresources(&limdfsresources,dfst);
    }
  }
  FREEARRAY(&storeonline,Simplematch);
  FREEARRAY(&storeoffline,Simplematch);
  if (mor != NULL)
  {
    freeMyersonlineresources(&mor);
  }
  if (genericindex == NULL)
  {
    gt_assert(encseq != NULL);
    freeEncodedsequence((Encodedsequence **) &encseq);
  } else
  {
    genericindex_delete(genericindex);
  }
  gt_seqiterator_delete(seqit);
  freeverboseinfo(&verboseinfo);
  return haserr ? -1 : 0;
}
