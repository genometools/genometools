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

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <errno.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "core/arraydef.h"
#include "core/assert_api.h"
#include "core/error_api.h"
#include "core/unused_api.h"
#include "core/progressbar.h"
#include "core/minmax.h"
#include "core/fa.h"
#include "core/progress_timer_api.h"
#include "core/encseq.h"
#include "core/safecast-gen.h"
#include "core/log_api.h"
#include "core/mathsupport.h"
#include "core/spacecalc.h"
#include "intcode-def.h"
#include "esa-fileend.h"
#include "sfx-diffcov.h"
#include "sfx-partssuf.h"
#include "sfx-suffixer.h"
#include "sfx-enumcodes.h"
#include "sfx-strategy.h"
#include "sfx-copysort.h"
#include "sfx-mappedstr.h"
#include "sfx-bentsedg.h"
#include "sfx-suffixgetset.h"
#include "stamp.h"

typedef struct
{
  unsigned long allocatedSuffixptr, nextfreeSuffixptr;
  GtSuffixsortspace *sssp;
} GtSuffixposbuffer;

struct Sfxiterator
{
  bool storespecials;
  GtCodetype currentmincode,
           currentmaxcode;
  unsigned long specialcharacters,
                widthofpart,
                totallength;
  GtSuffixsortspace *suffixsortspace;
  unsigned long nextfreeCodeatposition;
  Codeatposition *spaceCodeatposition;
  Suftabparts *suftabparts;
  const GtEncseq *encseq;
  GtReadmode readmode;
  Outlcpinfo *outlcpinfo;
  unsigned int part,
               numofchars,
               prefixlength;
  GtSuffixposbuffer fusp;
  GtSpecialrangeiterator *sri;
  GtRange overhang;
  bool exhausted;
  Bcktab *bcktab;
  GtCodetype numofallcodes;
  unsigned long *leftborder; /* points to bcktab->leftborder */
  unsigned long long bucketiterstep; /* for progressbar */
  Sfxstrategy sfxstrategy;
  GtLogger *logger;
  GtProgressTimer *sfxprogress;
  bool withprogressbar;
  Differencecover *dcov;
};

#ifdef SKDEBUG
static unsigned long iterproduceCodeatposition(Codeatposition *codelist,
                                               const  GtEncseq *encseq,
                                               GtReadmode readmode,
                                               unsigned int prefixlength,
                                               unsigned int numofchars)
{
  if (prefixlength > 1U)
  {
    Enumcodeatposition *ecp;
    unsigned long insertindex;
    Specialcontext specialcontext;

    ecp = gt_Enumcodeatposition_new(encseq,
                                    readmode,
                                    prefixlength,
                                    numofchars);
    for (insertindex = 0; gt_Enumcodeatposition_next(&specialcontext,ecp);
         insertindex++)
    {
      codelist[insertindex].maxprefixindex = specialcontext.maxprefixindex;
      codelist[insertindex].position = specialcontext.position;
      codelist[insertindex].code
        = gt_Enumcodeatposition_filledqgramcode(ecp,
                                                specialcontext.maxprefixindex,
                                                specialcontext.position -
                                                specialcontext.maxprefixindex);
    }
    gt_Enumcodeatposition_delete(ecp);
    ecp = NULL;
    return insertindex;
  }
  return 0;
}
#endif

static void compareCodeatpositionlists(const Codeatposition *codelist1,
                                       unsigned long len1,
                                       const Codeatposition *codelist2,
                                       unsigned long len2)
{
  unsigned long idx;

  if (len1 != len2)
  {
    fprintf(stderr,"len1 = %lu != %lu = len2\n",len1,len2);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  for (idx=0; idx<len1; idx++)
  {
    if (codelist1[idx].position != codelist2[idx].position)
    {
      fprintf(stderr,"listlength = %lu,idx %lu, codelist1.position = %lu != %lu"
                     " = codelist2.position\n",len1,idx,
                      codelist1[idx].position,
                      codelist2[idx].position);
      fprintf(stderr,"codelist1.maxprefixindex = %u,%u="
                     "codelist2.maxprefixindex\n",
                      codelist1[idx].maxprefixindex,
                      codelist2[idx].maxprefixindex);
      fprintf(stderr,"codelist1.code = %u,%u = "
                     "codelist2.code\n",
                      codelist1[idx].code,
                      codelist2[idx].code);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    if (codelist1[idx].maxprefixindex != codelist2[idx].maxprefixindex)
    {
      fprintf(stderr,"idx %lu, codelist1.maxprefixindex = %u != %u = "
                     "codelist2.maxprefixindex\n",idx,
                      codelist1[idx].maxprefixindex,
                      codelist2[idx].maxprefixindex);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    if (codelist1[idx].code != codelist2[idx].code)
    {
      fprintf(stderr,"idx %lu, codelist1.code = %u != %u = "
                     "codelist2.code\n",idx,
                      codelist1[idx].code,
                      codelist2[idx].code);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
}

#ifdef SKDEBUG
static void verifycodelistcomputation(
                       const GtEncseq *encseq,
                       GtReadmode readmode,
                       unsigned long realspecialranges,
                       unsigned int prefixlength,
                       unsigned int numofchars,
                       unsigned long nextfreeCodeatposition1,
                       const Codeatposition *spaceCodeatposition1)
{
  unsigned long nextfreeCodeatposition2;
  Codeatposition *spaceCodeatposition2;

  spaceCodeatposition2 = gt_malloc(sizeof (*spaceCodeatposition2) *
                                   (realspecialranges+1));
  nextfreeCodeatposition2 = iterproduceCodeatposition(spaceCodeatposition2,
                                                      encseq,
                                                      readmode,
                                                      prefixlength,
                                                      numofchars);
  gt_assert(realspecialranges+1 >= (unsigned long) nextfreeCodeatposition2);
  compareCodeatpositionlists(spaceCodeatposition1,
                             nextfreeCodeatposition1,
                             spaceCodeatposition2,
                             nextfreeCodeatposition2);
  gt_free(spaceCodeatposition2);
}
#endif

#ifdef SKDEBUG
static GtCodetype getencseqcode(const GtEncseq *encseq,
                                GtReadmode readmode,
                                unsigned long totallength,
                                const GtCodetype **multimappower,
                                unsigned int prefixlength,
                                unsigned long pos)
{
  GtCodetype code = 0;
  unsigned int idx;
  GtUchar cc;

  for (idx=0; idx<prefixlength; idx++)
  {
    gt_assert((unsigned long) (pos + idx) < totallength);
    cc = gt_encseq_get_encoded_char_nospecial(encseq,pos + idx,readmode);
    gt_assert(ISNOTSPECIAL(cc));
    code += multimappower[idx][cc];
  }
  return code;
}

static GtCodetype previouscode = 0;
static bool previouskmercodedefined = false,
            previousstorespecials = false;
unsigned int previousspecialpos = 0;
#endif

static void updatekmercount(void *processinfo,
                            unsigned long position,
                            const GtKmercode *kmercode)
{
  Sfxiterator *sfi = (Sfxiterator *) processinfo;

  if (kmercode->definedspecialposition)
  {
    if (sfi->storespecials)
    {
      if (kmercode->specialposition > 0)
      {
        if (sfi->sfxstrategy.storespecialcodes)
        {
          Codeatposition *cp;

          cp = sfi->spaceCodeatposition + sfi->nextfreeCodeatposition++;
          gt_assert(kmercode->code <= (GtCodetype) MAXCODEVALUE);
          cp->code = (unsigned int) kmercode->code;
          gt_assert(kmercode->specialposition
                    <= (unsigned int) MAXPREFIXLENGTH);
          cp->maxprefixindex = kmercode->specialposition;
          cp->position = position + kmercode->specialposition;
          /*
          printf("store(code=%u,maxprefixindex=%u,pos=%lu)\n",
                  cp->code,cp->maxprefixindex,cp->position);
          */
        }
        sfi->storespecials = false;
        gt_assert(kmercode->code > 0);
        sfi->leftborder[kmercode->code]++;
      }
    } else
    {
      if (kmercode->specialposition > 0)
      {
        gt_assert(kmercode->code > 0);
        sfi->leftborder[kmercode->code]++;
      } else
      {
        sfi->storespecials = true;
      }
    }
  } else
  {
#ifdef SKDEBUG
    if (kmercode->code == 0)
    {
      GtCodetype code2 = getencseqcode(sfi->encseq,
                                       sfi->readmode,
                                       gt_encseq_total_length(sfi->encseq),
                                       gt_bcktab_multimappower(sfi->bcktab),
                                       sfi->prefixlength,
                                       position);
      if (code2 != 0)
      {
        fprintf(stderr,"### position %lu, code2 = %lu != 0\n",position,code2);
        fprintf(stderr,"previouscode = " FormatGtCodetype "\n",previouscode);
        if (previouskmercodedefined)
        {
          fprintf(stderr,"previouskmercodedefined = true\n");
          fprintf(stderr,"previousstorespecials = %s\n",
                  previousstorespecials ? "true" : "false");
          fprintf(stderr,"previousspecialpos = %u\n",previousspecialpos);
        }
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
    }
#endif
    sfi->leftborder[kmercode->code]++;
  }
#ifdef SKDEBUG
  previouscode = kmercode->code;
  previouskmercodedefined = kmercode->defined;
  previousstorespecials = sfi->storespecials;
  previousspecialpos = kmercode->specialposition;
#endif
}

static void insertwithoutspecial(void *processinfo,
                                 unsigned long position,
                                 const GtKmercode *kmercode)
{
  if (!kmercode->definedspecialposition)
  {
    Sfxiterator *sfi = (Sfxiterator *) processinfo;

    if (kmercode->code >= sfi->currentmincode &&
        kmercode->code <= sfi->currentmaxcode)
    {
      unsigned long stidx = --sfi->leftborder[kmercode->code];
      gt_suffixsortspace_setdirectwithoffset(sfi->suffixsortspace,stidx,
                                             position);
      /* from right to left */
    }
  }
}

static void reversespecialcodes(Codeatposition *spaceCodeatposition,
                                unsigned long nextfreeCodeatposition)
{
  Codeatposition *front, *back, tmp;

  for (front = spaceCodeatposition,
       back = spaceCodeatposition + nextfreeCodeatposition - 1;
       front < back; front++, back--)
  {
    tmp = *front;
    *front = *back;
    *back = tmp;
  }
}

static void sfx_derivespecialcodesfromtable(Sfxiterator *sfi,bool deletevalues)
{
  GtCodetype code;
  unsigned int prefixindex;
  unsigned long insertindex, j, stidx;

  for (prefixindex=1U; prefixindex < sfi->prefixlength; prefixindex++)
  {
    for (j=0, insertindex = 0; j < sfi->nextfreeCodeatposition; j++)
    {
      if (prefixindex <= sfi->spaceCodeatposition[j].maxprefixindex)
      {
        code = gt_codedownscale(sfi->bcktab,
                                (GtCodetype) sfi->spaceCodeatposition[j].code,
                                prefixindex,
                                sfi->spaceCodeatposition[j].maxprefixindex);
        if (code >= sfi->currentmincode && code <= sfi->currentmaxcode)
        {
          gt_updatebckspecials(sfi->bcktab,code,sfi->numofchars,prefixindex);
          stidx = --sfi->leftborder[code];
          /* from right to left */
          gt_suffixsortspace_setdirectwithoffset(sfi->suffixsortspace,stidx,
                           sfi->spaceCodeatposition[j].position - prefixindex);
        }
      }
      if (deletevalues)
      {
        if (prefixindex < sfi->prefixlength - 1 &&
            prefixindex < sfi->spaceCodeatposition[j].maxprefixindex)
        {
          if (insertindex < j)
          {
            sfi->spaceCodeatposition[insertindex] =
              sfi->spaceCodeatposition[j];
          }
          insertindex++;
        }
      }
    }
    if (deletevalues)
    {
      sfi->nextfreeCodeatposition = insertindex;
    }
  }
}

static void sfx_derivespecialcodesonthefly(Sfxiterator *sfi)
{
  GtCodetype code;
  unsigned int prefixindex;
  unsigned long stidx;
  Enumcodeatposition *ecp;
  Specialcontext specialcontext;

  for (prefixindex=1U; prefixindex < sfi->prefixlength; prefixindex++)
  {
    ecp = gt_Enumcodeatposition_new(sfi->encseq,sfi->readmode,
                                    sfi->prefixlength,
                                    sfi->numofchars);
    while (gt_Enumcodeatposition_next(&specialcontext,ecp))
    {
      if (prefixindex <= specialcontext.maxprefixindex)
      {
        if (gt_Enumcodeatposition_filledqgramcodestopatmax(
                                            &code,
                                            ecp,
                                            prefixindex,
                                            specialcontext.position-prefixindex,
                                            sfi->currentmaxcode))
        {
          gt_assert(code <= sfi->currentmaxcode);
          if (code >= sfi->currentmincode)
          {
            gt_updatebckspecials(sfi->bcktab,code,sfi->numofchars,prefixindex);
            gt_assert(code > 0);
            stidx = --sfi->leftborder[code];
            /* from right to left */
            gt_suffixsortspace_setdirectwithoffset(sfi->suffixsortspace,stidx,
                          specialcontext.position - prefixindex);
          }
        }
      }
    }
    gt_Enumcodeatposition_delete(ecp);
    ecp = NULL;
  }
}

void gt_Sfxiterator_delete(Sfxiterator *sfi)
{
  if (sfi == NULL)
  {
    return;
  }
#ifdef SKDEBUG
  if (sfi->bcktab != NULL)
  {
    checkcountspecialcodes(sfi->bcktab);
  }
#endif
  if (sfi->bcktab != NULL)
  {
    gt_addfinalbckspecials(sfi->bcktab,sfi->numofchars,sfi->specialcharacters);
  }
  if (sfi->sri != NULL)
  {
    gt_specialrangeiterator_delete(sfi->sri);
  }
  gt_free(sfi->spaceCodeatposition);
  sfi->spaceCodeatposition = NULL;
  gt_suffixsortspace_delete(sfi->suffixsortspace,true);
  gt_freesuftabparts(sfi->suftabparts);
  gt_bcktab_delete(sfi->bcktab);
  if (sfi->dcov != NULL)
  {
    gt_differencecover_delete(sfi->dcov);
  }
  gt_free(sfi);
}

#ifdef SKDEBUG
static void showleftborder(const unsigned long *leftborder,
                           GtCodetype numofallcodes)
{
  GtCodetype i;

  for (i=0; i<MIN(numofallcodes,(GtCodetype) 1024); i++)
  {
    printf("leftborder[" FormatGtCodetype "]=%lu\n",
            i,leftborder[i]);
  }
}
#endif

static void getencseqkmersupdatekmercount(const GtEncseq *encseq,
                                          GtReadmode readmode,
                                          unsigned int kmersize,
                                          Sfxiterator *sfi)
{
  GtKmercodeiterator *kmercodeiterator;
  const GtKmercode *kmercodeptr;

  kmercodeiterator = gt_kmercodeiterator_encseq_new(encseq,readmode,kmersize,0);
  if (!gt_kmercodeiterator_inputexhausted(kmercodeiterator))
  {
    unsigned long position = 0;

    while ((kmercodeptr = gt_kmercodeiterator_encseq_next(kmercodeiterator))
                          != NULL)
    {
      updatekmercount(sfi,position++,kmercodeptr);
    }
  }
  gt_kmercodeiterator_delete(kmercodeiterator);
}

void getencseqkmersinsertwithoutspecial(const GtEncseq *encseq,
                                        GtReadmode readmode,
                                        unsigned int kmersize,
                                        Sfxiterator *sfi)
{
  GtKmercodeiterator *kmercodeiterator;
  const GtKmercode *kmercodeptr;

  kmercodeiterator = gt_kmercodeiterator_encseq_new(encseq,readmode,kmersize,0);
  if (!gt_kmercodeiterator_inputexhausted(kmercodeiterator))
  {
    unsigned long position = 0;

    while ((kmercodeptr = gt_kmercodeiterator_encseq_next(kmercodeiterator))
                          != NULL)
    {
      insertwithoutspecial(sfi,position++,kmercodeptr);
    }
    gt_kmercodeiterator_delete(kmercodeiterator);
  }
}

static int computepartsfittingmaximumspace(size_t estimatedspace,
                                           unsigned long maximumspace,
                                           const unsigned long *leftborder,
                                           unsigned long numofallcodes,
                                           unsigned long totallength,
                                           unsigned long specialcharacters,
                                           bool suftabasulongarray,
                                           GtError *err)
{
  unsigned int parts;
  Suftabparts *suftabparts;

  if (estimatedspace >= (size_t) maximumspace)
  {
    gt_error_set(err,"already used %.2f MB of memory cannot compute "
                     "enhanced suffix array in at most %lu MB",
                     GT_MEGABYTES(estimatedspace),maximumspace);
    return -1;
  }
  for (parts = 1U; parts <= 100U; parts++)
  {
    size_t suftabsize;

    suftabparts = gt_newsuftabparts(parts,
                                    leftborder,
                                    numofallcodes,
                                    totallength - specialcharacters,
                                    specialcharacters + 1,
                                    NULL);
    gt_assert(suftabparts != NULL);
    suftabsize = gt_suffixsortspace_requiredspace(
                                              stpgetlargestwidth(suftabparts),
                                              totallength,
                                              suftabasulongarray);
    if ((unsigned long) (suftabsize + estimatedspace) <= maximumspace)
    {
      gt_freesuftabparts(suftabparts);
      return (int) parts;
    }
    gt_freesuftabparts(suftabparts);
  }
  gt_error_set(err,"cannot compute enhanced suffix array in at most %lu MB",
                   maximumspace);
  return -1;
}

#ifdef DEBUGSIZEESTIMATION
static void verifyestimatedspace(size_t estimatedspace)
{
  unsigned long usedspace_ma_fa = gt_ma_get_space_current() +
                                  gt_fa_get_space_current();
  if (usedspace_ma_fa > 0)
  {
    double relativedifference;

    if (usedspace_ma_fa >= (unsigned long) estimatedspace)
    {
      relativedifference
        = (double) (usedspace_ma_fa - estimatedspace)/usedspace_ma_fa;
    } else
    {
      relativedifference
        = (double) (estimatedspace - usedspace_ma_fa)/estimatedspace;
    }
    if (usedspace_ma_fa > 100000UL &&
        gt_double_larger_double(relativedifference,0.1))
    {
      fprintf(stderr, "relativedifference %.2f too large: estimatedspace=%.2f, "
                      "usedspace_ma_fa=%.2f\n",relativedifference,
                                               GT_MEGABYTES(estimatedspace),
                                               GT_MEGABYTES(usedspace_ma_fa));
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
}
#endif

static void gt_sfxcheckkmercode(void *processinfo,
                                GT_UNUSED unsigned long pos,
                                GtCodetype code)
{
  GtKmercodeiterator *kmercodeiterator = (GtKmercodeiterator *) processinfo;
  const GtKmercode *kmercodeptr;

  gt_assert(kmercodeiterator != NULL);
  kmercodeptr = gt_kmercodeiterator_encseq_nonspecial_next(kmercodeiterator);
  gt_assert(kmercodeptr != NULL && !kmercodeptr->definedspecialposition);
  gt_assert(code == kmercodeptr->code);
}

typedef struct
{
  unsigned long nextfreeCodeatposition;
  Codeatposition *spaceCodeatposition;
} GtArrayCodeatposition;

static void gt_storespecialcode(void *info,
                                unsigned int maxprefixindex,
                                unsigned int code,
                                unsigned long position)
{
  GtArrayCodeatposition *arrspca = (GtArrayCodeatposition *) info;
  Codeatposition *spcaptr;

  spcaptr = arrspca->spaceCodeatposition + arrspca->nextfreeCodeatposition++;
  spcaptr->maxprefixindex = maxprefixindex;
  spcaptr->code = code;
  spcaptr->position = position;
}

Sfxiterator *gt_Sfxiterator_new(const GtEncseq *encseq,
                                GtReadmode readmode,
                                unsigned int prefixlength,
                                unsigned int numofparts,
                                unsigned long maximumspace,
                                void *voidoutlcpinfo,
                                const Sfxstrategy *sfxstrategy,
                                GtProgressTimer *sfxprogress,
                                bool withprogressbar,
                                GtLogger *logger,
                                GtError *err)
{
  Sfxiterator *sfi = NULL;
  unsigned long realspecialranges, specialcharacters;
  size_t estimatedspace;
  bool haserr = false;

  gt_error_check(err);
  /*
  printf("spacepeak at start of gt_Sfxiterator_new = %lu\n",
          gt_ma_get_space_peak());*/ /* in bytes */
  gt_assert(encseq != NULL);
  estimatedspace = (size_t) gt_encseq_sizeofrep(encseq);
  realspecialranges = gt_encseq_realspecialranges(encseq);
  specialcharacters = gt_encseq_specialcharacters(encseq);
  gt_assert(prefixlength > 0);
  if (sfxstrategy != NULL && sfxstrategy->storespecialcodes &&
      prefixlength > (unsigned int) MAXPREFIXLENGTH)
  {
    gt_error_set(err,"argument for option -pl must be in the range [1,%u]",
                  MAXPREFIXLENGTH);
    haserr = true;
  }
  if (!haserr)
  {
    sfi = gt_malloc(sizeof (*sfi));
    if (sfxstrategy != NULL && sfxstrategy->storespecialcodes)
    {
      sfi->spaceCodeatposition
        = gt_malloc(sizeof (*sfi->spaceCodeatposition) * (realspecialranges+1));
      gt_log_log("sizeof (spaceCodeatposition)=%lu bytes",
                 (unsigned long) (sizeof (*sfi->spaceCodeatposition) *
                                          (realspecialranges+1)));
      estimatedspace += sizeof (*sfi->spaceCodeatposition) *
                        (realspecialranges+1);
    } else
    {
      sfi->spaceCodeatposition = NULL;
    }
    sfi->bcktab = NULL;
    sfi->nextfreeCodeatposition = 0;
    sfi->suffixsortspace = NULL;
    sfi->suftabparts = NULL;
    sfi->encseq = encseq;
    sfi->readmode = readmode;
    sfi->numofchars = gt_alphabet_num_of_chars(gt_encseq_alphabet(encseq));
    sfi->prefixlength = prefixlength;
    sfi->dcov = NULL;
    sfi->withprogressbar = withprogressbar;
    if (sfxstrategy != NULL)
    {
       sfi->sfxstrategy = *sfxstrategy;
       if (sfxstrategy->cmpcharbychar || !gt_encseq_bitwise_cmp_ok(encseq))
       {
         sfi->sfxstrategy.cmpcharbychar = true;
       } else
       {
         sfi->sfxstrategy.cmpcharbychar = false;
       }
    } else
    {
      defaultsfxstrategy(&sfi->sfxstrategy,
                         gt_encseq_bitwise_cmp_ok(encseq) ? false : true);
    }
    gt_logger_log(logger,"maxinsertionsort=%lu",
                sfi->sfxstrategy.maxinsertionsort);
    gt_logger_log(logger,"maxbltriesort=%lu",
                sfi->sfxstrategy.maxbltriesort);
    gt_logger_log(logger,"maxcountingsort=%lu",
                sfi->sfxstrategy.maxcountingsort);
    gt_logger_log(logger,"storespecialcodes=%s",
                sfi->sfxstrategy.storespecialcodes ? "true" : "false");
    gt_logger_log(logger,"cmpcharbychar=%s",
                sfi->sfxstrategy.cmpcharbychar ? "true" : "false");
    if (sfi->sfxstrategy.ssortmaxdepth.defined)
    {
      gt_logger_log(logger,"ssortmaxdepth=%u",
                               sfi->sfxstrategy.ssortmaxdepth.valueunsignedint);
    } else
    {
      gt_logger_log(logger,"ssortmaxdepth=undefined");
    }
    sfi->totallength = gt_encseq_total_length(encseq);
    gt_logger_log(logger,"totallength=%lu",sfi->totallength);
    sfi->specialcharacters = specialcharacters;
    sfi->outlcpinfo = (Outlcpinfo *) voidoutlcpinfo;
    sfi->sri = NULL;
    sfi->part = 0;
    sfi->exhausted = false;
    sfi->bucketiterstep = 0;
    sfi->logger = logger;
    sfi->sfxprogress = sfxprogress;

    if (sfi->sfxstrategy.differencecover > 0 &&
        gt_encseq_specialcharacters(encseq)
          < gt_encseq_total_length(encseq))
    {
      if (sfxprogress != NULL)
      {
        gt_progress_timer_start_new_state(sfxprogress,
                                          "sorting difference cover sample",
                                          stdout);
      }
      sfi->dcov = gt_differencecover_new(sfi->sfxstrategy.differencecover,
                                         encseq,readmode,logger);
      if (sfi->dcov == NULL)
      {
        gt_error_set(err,"no difference cover modulo %u found",
                     sfi->sfxstrategy.differencecover);
        haserr = true;
      } else
      {
        if (gt_differencecover_vparamverify(sfi->dcov,err) != 0)
        {
          haserr = true;
          gt_differencecover_delete(sfi->dcov);
          sfi->dcov = NULL;
        } else
        {
          gt_logger_log(logger,"presorting sample suffixes according to "
                                  "difference cover modulo %u",
                                  sfi->sfxstrategy.differencecover);
          gt_differencecover_sortsample(sfi->dcov,
                                        sfi->sfxstrategy.cmpcharbychar,
                                        false);
          estimatedspace += gt_differencecover_requiredspace(sfi->dcov);
        }
      }
    }
  }
  if (!haserr)
  {
    gt_assert(sfi != NULL);
    sfi->bcktab = gt_allocBcktab(sfi->numofchars,
                                 prefixlength,
                                 sfi->sfxstrategy.storespecialcodes,
                                 err);
    estimatedspace
      += (size_t) gt_sizeofbuckettable(sfi->numofchars,prefixlength);
    if (sfi->bcktab == NULL)
    {
      haserr = true;
      sfi->leftborder = NULL;
      sfi->numofallcodes = 0;
    } else
    {
      sfi->leftborder = gt_bcktab_leftborder(sfi->bcktab);
      sfi->numofallcodes = gt_bcktab_numofallcodes(sfi->bcktab);
    }
  }
  if (!haserr)
  {
    unsigned long largestbucketsize;
    gt_assert(sfi != NULL);
    sfi->storespecials = true;
    if (sfxprogress != NULL)
    {
      gt_progress_timer_start_new_state(sfxprogress,
                                        "counting prefix distribution",
                                        stdout);
    }
    gt_assert(sfi->leftborder != NULL);
    if (prefixlength == 1U)
    {
      unsigned int charidx;

      for (charidx=0; charidx<sfi->numofchars; charidx++)
      {
        sfi->leftborder[GT_ISDIRCOMPLEMENT(readmode) ?
                        GT_COMPLEMENTBASE(charidx) :
                        charidx]
          = gt_encseq_charcount(encseq,(GtUchar) charidx);
      }
    } else
    {
      if (sfi->sfxstrategy.iteratorbasedkmerscanning)
      {
        getencseqkmersupdatekmercount(encseq, readmode, prefixlength, sfi);
      } else
      {
        getencseqkmers(encseq,readmode,prefixlength,updatekmercount,sfi);
      }
      if (sfi->sfxstrategy.storespecialcodes)
      {
        gt_assert(realspecialranges+1
                  >= (unsigned long) sfi->nextfreeCodeatposition);
        reversespecialcodes(sfi->spaceCodeatposition,
                            sfi->nextfreeCodeatposition);
      }
#ifdef SKDEBUG
      verifycodelistcomputation(encseq,
                                readmode,
                                realspecialranges,
                                prefixlength,
                                sfi->numofchars,
                                sfi->nextfreeCodeatposition,
                                sfi->spaceCodeatposition);
#endif
      if (gt_has_twobitencoding(encseq) && sfi->sfxstrategy.storespecialcodes)
      {
        GtArrayCodeatposition codelist;
        GtKmercodeiterator *kmercodeiterator
          = gt_kmercodeiterator_encseq_new(encseq,readmode,prefixlength,0);

        codelist.nextfreeCodeatposition = 0;
        codelist.spaceCodeatposition
          = gt_malloc(sizeof (*codelist.spaceCodeatposition) *
                      (realspecialranges+1));
        getencseqkmers_twobitencoding(encseq,
                                      readmode,
                                      prefixlength,
                                      gt_sfxcheckkmercode,
                                      kmercodeiterator,
                                      gt_storespecialcode,
                                      &codelist);
        gt_kmercodeiterator_delete(kmercodeiterator);
        gt_assert(sfi->spaceCodeatposition != NULL);
        reversespecialcodes(codelist.spaceCodeatposition,
                            codelist.nextfreeCodeatposition);
        compareCodeatpositionlists(sfi->spaceCodeatposition,
                                   sfi->nextfreeCodeatposition,
                                   codelist.spaceCodeatposition,
                                   codelist.nextfreeCodeatposition);
          /*printf("# compared codelist\n");*/
        gt_free(codelist.spaceCodeatposition);
      }
    }
#ifdef SKDEBUG
    showleftborder(sfi->leftborder,sfi->numofallcodes);
#endif
    largestbucketsize
      = gt_bcktab_leftborderpartialsums(sfi->bcktab,
                                        sfi->totallength - specialcharacters);
    estimatedspace += sizeof (uint8_t) * largestbucketsize;
#ifdef DEBUGSIZEESTIMATION
    verifyestimatedspace(estimatedspace);
#endif
    if (maximumspace > 0)
    {
      int retval;

      gt_assert(numofparts == 1U);
      retval
        = computepartsfittingmaximumspace(estimatedspace,
                                          maximumspace,
                                          sfi->leftborder,
                                          sfi->numofallcodes,
                                          sfi->totallength,
                                          specialcharacters,
                                          sfi->sfxstrategy.suftabasulongarray,
                                          err);
      if (retval < 0)
      {
        haserr = true;
      } else
      {
        gt_assert(retval > 0);
        numofparts = (unsigned int) retval;
        gt_logger_log(logger, "derived parts=%u",numofparts);
      }
    }
  }
  if (!haserr)
  {
    gt_assert(sfi != NULL);
    sfi->suftabparts = gt_newsuftabparts(numofparts,
                                         sfi->leftborder,
                                         sfi->numofallcodes,
                                         sfi->totallength - specialcharacters,
                                         specialcharacters + 1,
                                         logger);
    gt_assert(sfi->suftabparts != NULL);
    sfi->suffixsortspace
      = gt_suffixsortspace_new(stpgetlargestwidth(sfi->suftabparts),
                               sfi->totallength,
                               sfi->sfxstrategy.suftabasulongarray);
    if (gt_encseq_has_specialranges(sfi->encseq))
    {
      sfi->sri = gt_specialrangeiterator_new(sfi->encseq,
                                             GT_ISDIRREVERSE(sfi->readmode)
                                               ? false : true);
    } else
    {
      sfi->sri = NULL;
    }
    sfi->fusp.sssp = sfi->suffixsortspace;
    sfi->fusp.allocatedSuffixptr = stpgetlargestwidth(sfi->suftabparts);
    sfi->overhang.start = sfi->overhang.end = 0;
  }
  if (haserr)
  {
    gt_Sfxiterator_delete(sfi);
    return NULL;
  }
  return sfi;
}

static void preparethispart(Sfxiterator *sfi)
{
  unsigned long partwidth;
  unsigned int numofparts = stpgetnumofparts(sfi->suftabparts);

  if (sfi->part == 0 && sfi->withprogressbar)
  {
    gt_progressbar_start(&sfi->bucketiterstep,
                         (unsigned long long) sfi->numofallcodes);
  }
  sfi->currentmincode = stpgetcurrentmincode(sfi->part,sfi->suftabparts);
  sfi->currentmaxcode = stpgetcurrentmaxcode(sfi->part,sfi->suftabparts);
  sfi->widthofpart = stpgetcurrentwidthofpart(sfi->part,sfi->suftabparts);
  gt_suffixsortspace_offset_set(sfi->suffixsortspace,
                                stpgetcurrentsuftaboffset(sfi->part,
                                                          sfi->suftabparts));
  if (sfi->sfxstrategy.storespecialcodes)
  {
    sfx_derivespecialcodesfromtable(sfi,(numofparts == 1U) ? true : false);
  } else
  {
    sfx_derivespecialcodesonthefly(sfi);
  }
  if (sfi->sfxprogress != NULL)
  {
    gt_progress_timer_start_new_state(sfi->sfxprogress,
                                      "inserting suffixes into buckets",
                                      stdout);
  }
  if (sfi->sfxstrategy.iteratorbasedkmerscanning)
  {
    getencseqkmersinsertwithoutspecial(sfi->encseq,
                                       sfi->readmode,
                                       sfi->prefixlength,
                                       sfi);
  } else
  {
    getencseqkmers(sfi->encseq,sfi->readmode,sfi->prefixlength,
                   insertwithoutspecial,sfi);
  }
  if (sfi->sfxprogress != NULL)
  {
    gt_progress_timer_start_new_state(sfi->sfxprogress,
                                      "sorting the buckets",
                                      stdout);
  }
  /* exit(0); just for testing */
  partwidth = stpgetcurrentsumofwdith(sfi->part,sfi->suftabparts);
  if (sfi->sfxstrategy.ssortmaxdepth.defined &&
      sfi->prefixlength == sfi->sfxstrategy.ssortmaxdepth.valueunsignedint)
  {
    if (!sfi->sfxstrategy.streamsuftab)
    {
      /* option -maxdepth with argument */
      gt_qsufsort(sfi->suffixsortspace,
                  partwidth,
                  -1,
                  NULL,
                  sfi->encseq,
                  sfi->readmode,
                  sfi->currentmincode,
                  sfi->currentmaxcode,
                  sfi->bcktab,
                  sfi->numofchars,
                  sfi->prefixlength,
                  false,
                  true,
                  sfi->outlcpinfo);
    }
  } else
  {
    GtBucketspec2 *bucketspec2 = NULL;
    gt_assert(!sfi->sfxstrategy.streamsuftab);
    if (numofparts == 1U && sfi->outlcpinfo == NULL && sfi->prefixlength >= 2U)
    {
      bucketspec2 = gt_copysort_new(sfi->bcktab,sfi->encseq,sfi->readmode,
                                    partwidth,sfi->numofchars);
    }
    if (sfi->sfxstrategy.differencecover > 0)
    {
      gt_sortbucketofsuffixes(true,
                              sfi->suffixsortspace,
                              partwidth,
                              bucketspec2,
                              sfi->encseq,
                              sfi->readmode,
                              sfi->currentmincode,
                              sfi->currentmaxcode,
                              sfi->bcktab,
                              sfi->numofchars,
                              sfi->prefixlength,
                              &sfi->sfxstrategy,
                              (void *) sfi->dcov,
                              gt_differencecover_sortunsortedbucket,
                              sfi->logger);
    } else
    {
      gt_sortallbuckets (sfi->suffixsortspace,
                         bucketspec2,
                         sfi->encseq,
                         sfi->readmode,
                         sfi->currentmincode,
                         sfi->currentmaxcode,
                         partwidth,
                         sfi->bcktab,
                         sfi->numofchars,
                         sfi->prefixlength,
                         sfi->outlcpinfo,
                         &sfi->sfxstrategy,
                         &sfi->bucketiterstep,
                         sfi->logger);
    }
    if (bucketspec2 != NULL)
    {
      gt_copysort_derivesorting(bucketspec2,sfi->suffixsortspace,sfi->logger);
      gt_copysort_delete(bucketspec2);
      bucketspec2 = NULL;
    }
  }
  sfi->part++;
}

int gt_Sfxiterator_postsortfromstream(Sfxiterator *sfi, const GtStr *indexname,
                                      GtError *err)
{
  int mmapfiledesc = -1;
  GtStr *tmpfilename;
  struct stat sb;
  bool haserr = false;

  if (sfi->totallength == sfi->specialcharacters)
  {
    return 0;
  }
  gt_suffixsortspace_sortspace_delete(sfi->suffixsortspace);
  if (sfi->sfxstrategy.streamsuftab)
  {
    gt_assert(sfi->sfxstrategy.ssortmaxdepth.defined &&
              sfi->prefixlength ==
              sfi->sfxstrategy.ssortmaxdepth.valueunsignedint);
  }
  tmpfilename = gt_str_clone(indexname);
  gt_str_append_cstr(tmpfilename,SUFTABSUFFIX);
  mmapfiledesc = open(gt_str_get(tmpfilename), O_RDWR, 0);
  if (mmapfiledesc == -1)
  {
    gt_error_set(err,"cannot open file \"%s\": %s",gt_str_get(tmpfilename),
                 strerror(errno));
    haserr = true;
  }
  if (!haserr && fstat(mmapfiledesc, &sb) == -1)
  {
    gt_error_set(err,"cannot fstat file \"%s\": %s",gt_str_get(tmpfilename),
                 strerror(errno));
    haserr = true;
  }
  if (!haserr && sizeof (off_t) > sizeof (size_t) && sb.st_size > SIZE_MAX)
  {
    gt_error_set(err,"file \"%s\" of size %llu is too large to map",
                 gt_str_get(tmpfilename),(unsigned long long) sb.st_size);
    haserr = true;
  }
  if (!haserr
        && (size_t) sb.st_size != sizeof (unsigned long) * (sfi->totallength+1))
  {
    gt_error_set(err,"mapping file %s: file size "
                     " = %lu != %lu = expected number of units",
                 gt_str_get(tmpfilename),
                 (unsigned long) sb.st_size,
                 (unsigned long) (sfi->totallength+1) * sizeof (unsigned long));
    haserr = true;
  }
  if (!haserr)
  {
    gt_assert(sfi->totallength >= sfi->specialcharacters);
    gt_qsufsort(NULL,
                sfi->totallength - sfi->specialcharacters,
                mmapfiledesc,
                tmpfilename,
                sfi->encseq,
                sfi->readmode,
                0,
                sfi->currentmaxcode,
                sfi->bcktab,
                sfi->numofchars,
                sfi->prefixlength,
                sfi->sfxstrategy.hashexceptions,
                sfi->sfxstrategy.absoluteinversesuftab,
                sfi->outlcpinfo);
  }
  gt_str_delete(tmpfilename);
  if (close(mmapfiledesc) == -1)
  {
    gt_error_set(err,"cannot close file \"%s\": %s",gt_str_get(tmpfilename),
                 strerror(errno));
    haserr = true;
  }
  return haserr ? -1 : 0;
}

static void insertfullspecialrange(Sfxiterator *sfi,
                                   unsigned long leftpos,
                                   unsigned long rightpos)
{
  unsigned long pos;

  gt_assert(leftpos < rightpos);
  if (GT_ISDIRREVERSE(sfi->readmode))
  {
    pos = rightpos - 1;
  } else
  {
    pos = leftpos;
  }
  while (true)
  {
    if (GT_ISDIRREVERSE(sfi->readmode))
    {
      gt_assert(pos < sfi->totallength);
      gt_suffixsortspace_setdirect(sfi->fusp.sssp,sfi->fusp.nextfreeSuffixptr,
                                   GT_REVERSEPOS(sfi->totallength,pos));
      sfi->fusp.nextfreeSuffixptr++;
      if (pos == leftpos)
      {
        break;
      }
      pos--;
    } else
    {
      gt_suffixsortspace_setdirect(sfi->fusp.sssp,
                                   sfi->fusp.nextfreeSuffixptr,pos);
      sfi->fusp.nextfreeSuffixptr++;
      if (pos == rightpos-1)
      {
        break;
      }
      pos++;
    }
  }
}

static void fillspecialnextpage(Sfxiterator *sfi)
{
  GtRange range;
  unsigned long width;

  while (true)
  {
    if (sfi->overhang.start < sfi->overhang.end)
    {
      width = sfi->overhang.end - sfi->overhang.start;
      if (sfi->fusp.nextfreeSuffixptr + width > sfi->fusp.allocatedSuffixptr)
      {
        /* does not fit into the buffer, so only output a part */
        unsigned long rest = sfi->fusp.nextfreeSuffixptr +
                             width - sfi->fusp.allocatedSuffixptr;
        gt_assert(rest > 0);
        if (GT_ISDIRREVERSE(sfi->readmode))
        {
          insertfullspecialrange(sfi,sfi->overhang.start + rest,
                                 sfi->overhang.end);
          sfi->overhang.end = sfi->overhang.start + rest;
        } else
        {
          insertfullspecialrange(sfi,sfi->overhang.start,
                                     sfi->overhang.end - rest);
          sfi->overhang.start = sfi->overhang.end - rest;
        }
        break;
      }
      if (sfi->fusp.nextfreeSuffixptr + width == sfi->fusp.allocatedSuffixptr)
      { /* overhang fits into the buffer and buffer is full */
        insertfullspecialrange(sfi,sfi->overhang.start,
                               sfi->overhang.end);
        sfi->overhang.start = sfi->overhang.end = 0;
        break;
      }
      /* overhang fits into the buffer and buffer is not full */
      insertfullspecialrange(sfi,sfi->overhang.start,
                             sfi->overhang.end);
      sfi->overhang.start = sfi->overhang.end = 0;
    } else
    {
      if (sfi->sri != NULL && gt_specialrangeiterator_next(sfi->sri,&range))
      {
        width = range.end - range.start;
        gt_assert(width > 0);
        if (sfi->fusp.nextfreeSuffixptr + width > sfi->fusp.allocatedSuffixptr)
        { /* does not fit into the buffer, so only output a part */
          unsigned long rest = sfi->fusp.nextfreeSuffixptr +
                               width - sfi->fusp.allocatedSuffixptr;
          if (GT_ISDIRREVERSE(sfi->readmode))
          {
            insertfullspecialrange(sfi,range.start + rest,
                                   range.end);
            sfi->overhang.start = range.start;
            sfi->overhang.end = range.start + rest;
          } else
          {
            insertfullspecialrange(sfi,range.start,range.end - rest);
            sfi->overhang.start = range.end - rest;
            sfi->overhang.end = range.end;
          }
          break;
        }
        if (sfi->fusp.nextfreeSuffixptr + width == sfi->fusp.allocatedSuffixptr)
        { /* overhang fits into the buffer and buffer is full */
          insertfullspecialrange(sfi,range.start,range.end);
          sfi->overhang.start = sfi->overhang.end = 0;
          break;
        }
        insertfullspecialrange(sfi,range.start,range.end);
        sfi->overhang.start = sfi->overhang.end = 0;
      } else
      {
        if (sfi->fusp.nextfreeSuffixptr < sfi->fusp.allocatedSuffixptr)
        {
          gt_suffixsortspace_setdirect(sfi->fusp.sssp,
                                       sfi->fusp.nextfreeSuffixptr,
                                       sfi->totallength);
          sfi->fusp.nextfreeSuffixptr++;
          sfi->exhausted = true;
        }
        break;
      }
    }
  }
}

const GtSuffixsortspace *gt_Sfxiterator_next(unsigned long *numberofsuffixes,
                                             bool *specialsuffixes,
                                             Sfxiterator *sfi)
{
  if (sfi->part < stpgetnumofparts(sfi->suftabparts))
  {
    if (stpgetnumofparts(sfi->suftabparts) > 1U)
    {
      gt_logger_log(sfi->logger,"compute part %u",sfi->part+1);
    }
    preparethispart(sfi);
    *numberofsuffixes = sfi->widthofpart;
    *specialsuffixes = false;
    return sfi->suffixsortspace;
  }
  if (sfi->exhausted)
  {
    if (sfi->withprogressbar)
    {
      gt_progressbar_stop();
    }
    return NULL;
  }
  gt_suffixsortspace_offset_set(sfi->suffixsortspace,
                                sfi->totallength - sfi->specialcharacters);
  sfi->fusp.nextfreeSuffixptr = 0;
  fillspecialnextpage(sfi);
  gt_assert(sfi->fusp.nextfreeSuffixptr > 0);
  *numberofsuffixes = sfi->fusp.nextfreeSuffixptr;
  *specialsuffixes = true;
  return sfi->suffixsortspace;
}

int gt_Sfxiterator_bcktab2file(FILE *fp, const Sfxiterator *sfi, GtError *err)
{
  gt_error_check(err);
  return gt_bcktab2file(fp,sfi->bcktab,err);
}

unsigned long gt_Sfxiterator_longest(const Sfxiterator *sfi)
{
  return gt_suffixsortspace_longest(sfi->suffixsortspace);
}
