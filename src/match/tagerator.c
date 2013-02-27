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
#include "core/alphabet.h"
#include "core/arraydef.h"
#include "core/error.h"
#include "core/fileutils_api.h"
#include "core/format64.h"
#include "core/intbits.h"
#include "core/ma_api.h"
#include "core/seq_iterator_sequence_buffer_api.h"
#include "core/str_array.h"
#include "core/unused_api.h"
#include "apmeoveridx.h"
#include "dist-short.h"
#include "echoseq.h"
#include "esa-map.h"
#include "idx-limdfs.h"
#include "mssufpat.h"
#include "myersapm.h"
#include "revcompl.h"
#include "sarr-def.h"
#include "stamp.h"
#include "tagerator.h"

#define MAXTAGSIZE GT_INTWORDSIZE

#define ISRCDIR(TWL)  (((TWL)->tagptr == (TWL)->transformedtag)\
                        ? false\
                        : true)

typedef struct
{
  unsigned long dbstartpos,
         matchlength;
  bool rcmatch;
} TgrSimplematch;

typedef struct
{
  const GtUchar *tagptr;
  GtUchar transformedtag[MAXTAGSIZE],
        rctransformedtag[MAXTAGSIZE];
  unsigned long taglen;
} TgrTagwithlength;

typedef struct
{
  const TageratorOptions *tageratoroptions;
  unsigned int alphasize;
  const GtUchar *tagptr;
  const GtAlphabet *alpha;
  unsigned long *eqsvector;
  const TgrTagwithlength *twlptr;
  const GtEncseq *encseq;
} TgrShowmatchinfo;

#define ADDTABULATOR\
        if (firstitem)\
        {\
          firstitem = false;\
        } else\
        {\
          (void) putchar('\t');\
        }

static void tgr_showmatch(void *processinfo,const GtIdxMatch *match)
{
  TgrShowmatchinfo *showmatchinfo = (TgrShowmatchinfo *) processinfo;
  bool firstitem = true;

  gt_assert(showmatchinfo->tageratoroptions != NULL);
  if (showmatchinfo->tageratoroptions->outputmode & TAGOUT_DBLENGTH)
  {
    printf("%lu",match->dblen);
    firstitem = false;
  }
  if (showmatchinfo->tageratoroptions->outputmode & TAGOUT_DBSTARTPOS)
  {
    ADDTABULATOR;
    if (showmatchinfo->tageratoroptions->outputmode & TAGOUT_DBABSPOS)
    {
      printf("%lu",match->dbstartpos);
    } else
    {
      unsigned long seqstartpos,
                    seqnum = gt_encseq_seqnum(showmatchinfo->encseq,
                                                  match->dbstartpos);
      seqstartpos = gt_encseq_seqstartpos(showmatchinfo->encseq, seqnum);
      gt_assert(seqstartpos <= match->dbstartpos);
      printf("%lu\t%lu",seqnum, match->dbstartpos - seqstartpos);
    }
  }
  if (showmatchinfo->tageratoroptions->outputmode & TAGOUT_DBSEQUENCE)
  {
    ADDTABULATOR;
    gt_assert(match->dbsubstring != NULL);
    gt_alphabet_decode_seq_to_fp(showmatchinfo->alpha,
                                 stdout,
                                 match->dbsubstring,
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
          = gt_reversesuffixmatch(showmatchinfo->eqsvector,
                               showmatchinfo->alphasize,
                               match->dbsubstring,
                               (unsigned long) match->dblen,
                               showmatchinfo->tagptr,
                               match->querylen,
                               (unsigned long) showmatchinfo->tageratoroptions->
                                                        userdefinedmaxdistance);
        gt_assert(match->querylen >= suffixlength);
        if (showmatchinfo->tageratoroptions->outputmode & TAGOUT_TAGSTARTPOS)
        {
          ADDTABULATOR;
          printf("%lu",match->querylen - suffixlength);
        }
        if (showmatchinfo->tageratoroptions->outputmode & TAGOUT_TAGLENGTH)
        {
          ADDTABULATOR;
          printf("%lu",suffixlength);
        }
        if (showmatchinfo->tageratoroptions->outputmode & TAGOUT_TAGSUFFIXSEQ)
        {
          ADDTABULATOR;
          gt_alphabet_decode_seq_to_fp(NULL,
                                       stdout,showmatchinfo->tagptr +
                                       (match->querylen - suffixlength),
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
        printf("%lu",match->querylen);
      }
      if (showmatchinfo->tageratoroptions->outputmode & TAGOUT_TAGSUFFIXSEQ)
      {
        ADDTABULATOR;
        gt_alphabet_decode_seq_to_fp(NULL,
                                     stdout,
                                     showmatchinfo->tagptr,
                                     match->querylen);
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
  TgrSimplematch *spaceTgrSimplematch;
  unsigned long nextfreeTgrSimplematch, allocatedTgrSimplematch;
  const TgrTagwithlength *twlptr;
} ArrayTgrSimplematch;

static void tgr_storematch(void *processinfo,const GtIdxMatch *match)
{
  ArrayTgrSimplematch *storetab = (ArrayTgrSimplematch *) processinfo;
  TgrSimplematch *simplematch;

  GT_GETNEXTFREEINARRAY(simplematch,storetab,TgrSimplematch,32);
  simplematch->dbstartpos = match->dbstartpos;
  simplematch->matchlength = match->dblen;
  simplematch->rcmatch = ISRCDIR(storetab->twlptr);
}

static void checkmstats(void *processinfo,
                        const void *patterninfo,
                        unsigned long patternstartpos,
                        unsigned long mstatlength,
                        unsigned long leftbound,
                        unsigned long rightbound)
{
  GT_UNUSED unsigned long realmstatlength;
  TgrTagwithlength *twl = (TgrTagwithlength *) patterninfo;

  realmstatlength = genericmstats((const Limdfsresources *) processinfo,
                                  twl->tagptr + patternstartpos,
                                  twl->tagptr + twl->taglen);
#ifndef NDEBUG
  if (mstatlength != realmstatlength)
  {
    fprintf(stderr,"patternstartpos = %lu: mstatlength = %lu != %lu "
                   " = realmstatlength\n",
                    patternstartpos,mstatlength,realmstatlength);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
#endif
  if (gt_intervalwidthleq((const Limdfsresources *) processinfo,leftbound,
                       rightbound))
  {
    GtUchar cc;
    unsigned long *sptr, witnessposition;
    unsigned long idx;
    GtArrayGtUlong *mstatspos = gt_fromitv2sortedmatchpositions(
                                  (Limdfsresources *) processinfo,
                                  leftbound,
                                  rightbound,
                                  mstatlength);
    for (sptr = mstatspos->spaceGtUlong; sptr < mstatspos->spaceGtUlong +
                                               mstatspos->nextfreeGtUlong;
         sptr++)
    {
      witnessposition = *sptr;
      for (idx = patternstartpos; idx < patternstartpos + mstatlength; idx++)
      {
        cc = gt_limdfs_getencodedchar((const Limdfsresources *) processinfo,
                                  witnessposition + idx - patternstartpos,
                                  GT_READMODE_FORWARD);
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
          exit(GT_EXIT_PROGRAMMING_ERROR);
        }
      }
    }
  }
}

static void showmstats(void *processinfo,
                       const void *patterninfo,
                       GT_UNUSED unsigned long patternstartpos,
                       unsigned long mstatlength,
                       unsigned long leftbound,
                       unsigned long rightbound)
{
  TgrTagwithlength *twl = (TgrTagwithlength *) patterninfo;

  printf("%lu %c",mstatlength,ISRCDIR(twl) ? '-' : '+');
  if (gt_intervalwidthleq((const Limdfsresources *) processinfo,leftbound,
                       rightbound))
  {
    unsigned long idx;
    GtArrayGtUlong *mstatspos = gt_fromitv2sortedmatchpositions(
                                  (Limdfsresources *) processinfo,
                                  leftbound,
                                  rightbound,
                                  mstatlength);
    for (idx = 0; idx<mstatspos->nextfreeGtUlong; idx++)
    {
      printf(" %lu",mstatspos->spaceGtUlong[idx]);
    }
  }
  printf("\n");
}

static int cmpdescend(const void *a,const void *b)
{
  TgrSimplematch *valuea = (TgrSimplematch *) a;
  TgrSimplematch *valueb = (TgrSimplematch *) b;

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

static int dotransformtag(GtUchar *transformedtag,
                          const GtUchar *symbolmap,
                          const GtUchar *currenttag,
                          unsigned long taglen,
                          uint64_t tagnumber,
                          bool replacewildcard,
                          GtError *err)
{
  unsigned long idx;
  GtUchar charcode;

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
    if (charcode == (GtUchar) UNDEFCHAR)
    {
      gt_error_set(err,"undefined character '%c' in tag number " Formatuint64_t,
                currenttag[idx],
                PRINTuint64_tcast(tagnumber));
      return -1;
    }
    if (charcode == (GtUchar) WILDCARD)
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
                                 bool doonline,
                                 bool docompare,
                                 unsigned long maxintervalwidth,
                                 bool skpp,
                                 Myersonlineresources *mor,
                                 Limdfsresources *limdfsresources,
                                 const GtUchar *tagptr,
                                 unsigned long taglen)
{
  if (doonline || (!domstats && docompare))
  {
    gt_assert(mor != NULL);
    gt_edistmyersbitvectorAPM(mor,tagptr,taglen,maxdistance);
  }
  if (!doonline || docompare)
  {
    if (domstats)
    {
      gt_indexbasedmstats(limdfsresources,tagptr,taglen,dfst);
      return false;
    }
    if (maxdistance == 0)
    {
      return gt_indexbasedexactpatternmatching(limdfsresources,tagptr,taglen);
    } else
    {
      return gt_indexbasedapproxpatternmatching(limdfsresources,tagptr,taglen,
                                             maxdistance,
                                             maxintervalwidth,
                                             skpp,
                                             dfst);
    }
  }
  return false;
}

static void compareresults(const ArrayTgrSimplematch *storeonline,
                           const ArrayTgrSimplematch *storeoffline)
{
  unsigned long ss;

  if (storeonline->nextfreeTgrSimplematch
        != storeoffline->nextfreeTgrSimplematch)
  {
    fprintf(stderr,"nextfreeTgrSimplematch: storeonline = %lu != %lu "
                   "storeoffline\n",
                   storeonline->nextfreeTgrSimplematch,
                   storeoffline->nextfreeTgrSimplematch);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  gt_assert(storeonline->nextfreeTgrSimplematch ==
            storeoffline->nextfreeTgrSimplematch);
  if (storeoffline->nextfreeTgrSimplematch > 1UL)
  {
    qsort(storeoffline->spaceTgrSimplematch,(size_t)
          storeoffline->nextfreeTgrSimplematch,
          sizeof (TgrSimplematch),
          cmpdescend);
  }
  for (ss=0; ss < storeoffline->nextfreeTgrSimplematch; ss++)
  {
    gt_assert(storeonline->spaceTgrSimplematch != NULL &&
           storeoffline->spaceTgrSimplematch != NULL);
    if (storeonline->spaceTgrSimplematch[ss].rcmatch &&
        !storeoffline->spaceTgrSimplematch[ss].rcmatch)
    {
      fprintf(stderr,"rcmatch: storeonline[%lu] = p != d "
                     "= storeoffline[%lu]\n",ss,ss);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    if (!storeonline->spaceTgrSimplematch[ss].rcmatch &&
        storeoffline->spaceTgrSimplematch[ss].rcmatch)
    {
      fprintf(stderr,"rcmatch: storeonline[%lu] = d != p "
                     "= storeoffline[%lu]\n",ss,ss);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    if (storeonline->spaceTgrSimplematch[ss].matchlength !=
        storeoffline->spaceTgrSimplematch[ss].matchlength)
    {
      fprintf(stderr,"matchlength: storeonline[%lu] = %lu"
                     " != %lu = storeoffline[%lu]\n",
                     ss,
                     storeonline->spaceTgrSimplematch[ss].matchlength,
                     storeoffline->spaceTgrSimplematch[ss].matchlength,
                     ss);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    if (storeonline->spaceTgrSimplematch[ss].dbstartpos !=
        storeoffline->spaceTgrSimplematch[ss].dbstartpos)
    {
      fprintf(stderr,"dbstartpos: storeonline[%lu] = %lu"
                     " != %lu = storeoffline[%lu]\n",
                     ss,
                     storeonline->spaceTgrSimplematch[ss].dbstartpos,
                     storeoffline->spaceTgrSimplematch[ss].dbstartpos,
                     ss);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
}

static void searchoverstrands(const TageratorOptions *tageratoroptions,
                              TgrTagwithlength *twl,
                              const AbstractDfstransformer *dfst,
                              Myersonlineresources *mor,
                              Limdfsresources *limdfsresources,
                              TgrShowmatchinfo *showmatchinfo,
                              ArrayTgrSimplematch *storeonline,
                              ArrayTgrSimplematch *storeoffline)
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
                                 tageratoroptions->doonline,
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

int gt_runtagerator(const TageratorOptions *tageratoroptions,GtError *err)
{
  bool haserr = false, firstitem;
  int retval;
  Myersonlineresources *mor = NULL;
  Genericindex *genericindex = NULL;
  const GtEncseq *encseq = NULL;
  GtLogger *logger;

  logger = gt_logger_new(tageratoroptions->verbose,
                         GT_LOGGER_DEFLT_PREFIX, stdout);
  if (tageratoroptions->doonline)
  {
    GtEncseqLoader *el;
    el = gt_encseq_loader_new();
    gt_encseq_loader_do_not_require_des_tab(el);
    gt_encseq_loader_do_not_require_ssp_tab(el);
    gt_encseq_loader_do_not_require_sds_tab(el);
    gt_encseq_loader_set_logger(el, logger);
    encseq = gt_encseq_loader_load(el, gt_str_get(tageratoroptions->indexname),
                                   err);
    gt_encseq_loader_delete(el);
    if (encseq == NULL)
    {
      haserr = true;
    }
  } else
  {
    genericindex = genericindex_new(gt_str_get(tageratoroptions->indexname),
                                    tageratoroptions->withesa,
                                    tageratoroptions->withesa ||
                                    tageratoroptions->docompare,
                                    false,
                                    (tageratoroptions->outputmode &
                                     TAGOUT_DBABSPOS) ? false : true,
                                    tageratoroptions->userdefinedmaxdepth,
                                    logger,
                                    err);
    if (genericindex == NULL)
    {
      haserr = true;
    } else
    {
      encseq = genericindex_getencseq(genericindex);
    }
  }
  if (!haserr)
  {
    TgrTagwithlength twl;
    uint64_t tagnumber;
    unsigned int numofchars;
    const GtUchar *symbolmap, *currenttag;
    char *desc = NULL;
    ProcessIdxMatch processmatch;
    TgrShowmatchinfo showmatchinfo;
    void *processmatchinfoonline, *processmatchinfooffline;
    Limdfsresources *limdfsresources = NULL;
    const GtAlphabet *alpha;
    ArrayTgrSimplematch storeonline, storeoffline;
    const AbstractDfstransformer *dfst;
    GtSeqIterator *seqit = NULL;

    if (tageratoroptions->userdefinedmaxdistance >= 0)
    {
      dfst = gt_apme_AbstractDfstransformer();
    } else
    {
      dfst = gt_pms_AbstractDfstransformer();
    }
    GT_INITARRAY(&storeonline,TgrSimplematch);
    GT_INITARRAY(&storeoffline,TgrSimplematch);
    storeonline.twlptr = storeoffline.twlptr = &twl;
    alpha = gt_encseq_alphabet(encseq);
    symbolmap = gt_alphabet_symbolmap(alpha);
    numofchars = gt_alphabet_num_of_chars(alpha);
    if (tageratoroptions->docompare)
    {
      processmatch = tgr_storematch;
      processmatchinfoonline = &storeonline;
      processmatchinfooffline = &storeoffline;
      showmatchinfo.eqsvector = NULL;
      showmatchinfo.encseq = encseq;
    } else
    {
      processmatch = tgr_showmatch;
      showmatchinfo.twlptr = &twl;
      showmatchinfo.tageratoroptions = tageratoroptions;
      showmatchinfo.alphasize = (unsigned int) numofchars;
      showmatchinfo.alpha = alpha;
      showmatchinfo.eqsvector = gt_malloc(sizeof (*showmatchinfo.eqsvector) *
                                          showmatchinfo.alphasize);
      showmatchinfo.encseq = encseq;
      processmatchinfooffline = &showmatchinfo;
      processmatchinfoonline = &showmatchinfo;
    }
    if (tageratoroptions->doonline || tageratoroptions->docompare)
    {
      gt_assert(encseq != NULL);
      mor = gt_newMyersonlineresources(numofchars,
                                    tageratoroptions->nowildcards,
                                    encseq,
                                    processmatch,
                                    processmatchinfoonline);
    }
    if (!tageratoroptions->doonline || tageratoroptions->docompare)
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
      limdfsresources = gt_newLimdfsresources(genericindex,
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
    gt_getsetargmodekeywords(tageratoroptions->modedesc,
                             tageratoroptions->numberofmodedescentries,
                             tageratoroptions->outputmode);
    seqit = gt_seq_iterator_sequence_buffer_new(tageratoroptions->tagfiles,
                                                err);
    if (!seqit)
    {
      haserr = true;
    }
    if (!haserr)
    {
      for (tagnumber = 0; !haserr; tagnumber++)
      {
        retval = gt_seq_iterator_next(seqit, &currenttag, &twl.taglen, &desc,
                                     err);
        if (retval != 1)
        {
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
        gt_copy_reversecomplement(twl.rctransformedtag,twl.transformedtag,
                               twl.taglen);
        twl.tagptr = twl.transformedtag;
        firstitem = true;
        printf("#");
        if (tageratoroptions->outputmode & TAGOUT_TAGNUM)
        {
          printf("\t" Formatuint64_t,PRINTuint64_tcast(tagnumber));
          firstitem = false;
        }
        if (tageratoroptions->outputmode & TAGOUT_TAGLENGTH)
        {
          ADDTABULATOR;
          printf("%lu",twl.taglen);
        }
        if (tageratoroptions->outputmode & TAGOUT_TAGSEQ)
        {
          ADDTABULATOR;
          gt_alphabet_decode_seq_to_fp(alpha,stdout,twl.transformedtag,
                                       twl.taglen);
        }
        printf("\n");
        storeoffline.nextfreeTgrSimplematch = 0;
        storeonline.nextfreeTgrSimplematch = 0;
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
      }
      gt_seq_iterator_delete(seqit);
    }
    GT_FREEARRAY(&storeonline,TgrSimplematch);
    GT_FREEARRAY(&storeoffline,TgrSimplematch);
    gt_free(showmatchinfo.eqsvector);
    if (limdfsresources != NULL)
    {
      gt_freeLimdfsresources(&limdfsresources,dfst);
    }
  }
  gt_freeMyersonlineresources(mor);
  if (genericindex == NULL)
  {
    if (encseq != NULL)
    {
      gt_encseq_delete((GtEncseq *) encseq);
      encseq = NULL;
    }
  } else
  {
    genericindex_delete(genericindex);
  }
  gt_logger_delete(logger);
  return haserr ? -1 : 0;
}
