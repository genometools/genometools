/*
  Copyright (c) 2007-2013 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2013 Center for Bioinformatics, University of Hamburg

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
#ifndef S_SPLINT_S
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#endif
#include "core/arraydef.h"
#include "core/assert_api.h"
#include "core/chardef.h"
#include "core/error_api.h"
#include "core/unused_api.h"
#include "core/progressbar.h"
#include "core/minmax.h"
#include "core/fa.h"
#include "core/timer_api.h"
#include "core/encseq.h"
#include "core/safecast-gen.h"
#include "core/log_api.h"
#include "core/mathsupport.h"
#include "core/spacecalc.h"
#include "core/divmodmul.h"
#include "core/format64.h"
#include "core/fileutils.h"
#ifdef GT_THREADS_ENABLED
#include "core/thread_api.h"
#endif
#include "intcode-def.h"
#include "firstcodes-buf.h"
#include "esa-fileend.h"
#include "kmercodes.h"
#include "sfx-diffcov.h"
#include "sfx-partssuf.h"
#include "sfx-suffixer.h"
#include "sfx-enumcodes.h"
#include "sfx-strategy.h"
#include "sfx-copysort.h"
#include "sfx-mappedstr.h"
#include "sfx-bentsedg.h"
#include "sfx-suffixgetset.h"
#include "sfx-maprange.h"

struct Sfxiterator
{
  /* globally constant */
  const GtEncseq *encseq;
  GtReadmode readmode;
  GtUword specialcharacters,
                totallength;
  unsigned int numofchars,
               prefixlength;
  Sfxstrategy sfxstrategy;
  bool withprogressbar;

  /* invariant for each part */
  GtSuftabparts *suftabparts;
  GtOutlcpinfo *outlcpinfoforsample;
  GtBcktab *bcktab;
  GtLeftborder *leftborder; /* points to bcktab->leftborder */
  GtDifferencecover *dcov;

  /* changed in each part */
  GtSuffixsortspace *suffixsortspace;
  GtCodetype currentmincode,
             currentmaxcode;
  GtUword widthofpart;
  unsigned int part;
  GtOutlcpinfo *outlcpinfo;
  bool exhausted;
  GtUint64 bucketiterstep; /* for progressbar */
  GtLogger *logger;
  GtTimer *sfxprogress;
  GtSSSPbuf *sssp_buf;
  GtSpecialrangeiterator *sri; /* refers to space used in each part */

  /* use for generating k-mer codes */
  FILE *outfpbcktab;
  bool storespecials;
  GtUword nextfreeCodeatposition;
  Codeatposition *spaceCodeatposition;
  unsigned int kmerfastmaskright;
  GtSuffixsortspace_exportptr *exportptr;
  unsigned int spmopt_kmerscansize,
               spmopt_kmerscancodeshift2bckcode,
               spmopt_kmerscancodeshift2prefixcode,
               spmopt_additionalprefixchars;
  GtBitsequence *markprefixbuckets,
                *marksuffixbuckets;
  GtCodetype spmopt_kmerscancodesuffixmask;
  GtUword spmopt_numofallprefixcodes,
                spmopt_numofallsuffixcodes;
  GtSfxmappedrange *mappedmarkprefixbuckets;
#ifdef GT_THREADS_ENABLED
#ifdef GT_THREADS_PARTITION
  GtSuftabparts **partitions_for_threads;
#endif
#endif
};

#ifdef SKDEBUG
static GtUword iterproduceCodeatposition(Codeatposition *codelist,
                                               const  GtEncseq *encseq,
                                               GtReadmode readmode,
                                               unsigned int prefixlength,
                                               unsigned int numofchars)
{
  if (prefixlength > 1U)
  {
    Enumcodeatposition *ecp;
    GtUword insertindex;
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

static void compareCodeatpositionlists(const Codeatposition *codelist1,
                                       GtUword len1,
                                       const Codeatposition *codelist2,
                                       GtUword len2)
{
  GtUword idx;

  if (len1 != len2)
  {
    fprintf(stderr,"%s: len1 = "GT_WU" != "GT_WU" = len2\n",__func__,len1,len2);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  for (idx=0; idx<len1; idx++)
  {
    if (codelist1[idx].position != codelist2[idx].position)
    {
      fprintf(stderr,"%s: listlength = "GT_WU",idx "GT_WU", codelist1.position "
                     "= "GT_WU" != "GT_WU" = codelist2.position\n",
                      __func__,len1,idx,
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
      fprintf(stderr,"%s: idx "GT_WU", codelist1.maxprefixindex = %u != %u = "
                     "codelist2.maxprefixindex\n",__func__,idx,
                      codelist1[idx].maxprefixindex,
                      codelist2[idx].maxprefixindex);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    if (codelist1[idx].code != codelist2[idx].code)
    {
      fprintf(stderr,"%s: idx "GT_WU", codelist1.code = %u != %u = "
                     "codelist2.code\n",__func__,idx,
                      codelist1[idx].code,
                      codelist2[idx].code);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
}

static void verifycodelistcomputation(
                       const GtEncseq *encseq,
                       GtReadmode readmode,
                       GtUword realspecialranges,
                       unsigned int prefixlength,
                       unsigned int numofchars,
                       GtUword nextfreeCodeatposition1,
                       const Codeatposition *spaceCodeatposition1)
{
  GtUword nextfreeCodeatposition2;
  Codeatposition *spaceCodeatposition2;

  spaceCodeatposition2 = gt_malloc(sizeof (*spaceCodeatposition2) *
                                   (realspecialranges+1));
  nextfreeCodeatposition2 = iterproduceCodeatposition(spaceCodeatposition2,
                                                      encseq,
                                                      readmode,
                                                      prefixlength,
                                                      numofchars);
  gt_assert(realspecialranges+1 >= nextfreeCodeatposition2);
  compareCodeatpositionlists(spaceCodeatposition1,
                             nextfreeCodeatposition1,
                             spaceCodeatposition2,
                             nextfreeCodeatposition2);
  gt_free(spaceCodeatposition2);
}

static GtCodetype getencseqcode(const GtEncseq *encseq,
                                GtReadmode readmode,
                                GtUword totallength,
                                const GtCodetype **multimappower,
                                unsigned int prefixlength,
                                GtUword pos)
{
  GtCodetype code = 0;
  unsigned int idx;
  GtUchar cc;

  for (idx=0; idx<prefixlength; idx++)
  {
    gt_assert((GtUword) (pos + idx) < totallength);
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
                            GtUword position,
                            const GtKmercode *kmercode)
{
  Sfxiterator *sfi = (Sfxiterator *) processinfo;

  gt_assert(sfi->sfxstrategy.spmopt_minlength == 0);
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
          gt_assert(kmercode->code <= (GtCodetype) GT_MAXCODEVALUE);
          cp->code = (unsigned int) kmercode->code;
          gt_assert(kmercode->specialposition
                    <= (unsigned int) GT_MAXPREFIXLENGTH);
          cp->maxprefixindex = kmercode->specialposition;
          cp->position = position + kmercode->specialposition;
          /*
          printf("store(code=%u,maxprefixindex=%u,pos="GT_WU")\n",
                  cp->code,cp->maxprefixindex,cp->position);
          */
        }
        sfi->storespecials = false;
        gt_assert(kmercode->code > 0);
        gt_bcktab_leftborder_addcode(sfi->leftborder,kmercode->code);
      }
    } else
    {
      if (kmercode->specialposition > 0)
      {
        gt_assert(kmercode->code > 0);
        gt_bcktab_leftborder_addcode(sfi->leftborder,kmercode->code);
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
        fprintf(stderr,"%s: ### position "GT_WU", code2 = "GT_WU" != 0\n",
                __func__, position,code2);
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
    gt_bcktab_leftborder_addcode(sfi->leftborder,kmercode->code);
  }
#ifdef SKDEBUG
  previouscode = kmercode->code;
  previouskmercodedefined = kmercode->defined;
  previousstorespecials = sfi->storespecials;
  previousspecialpos = kmercode->specialposition;
#endif
}

#define GT_SCANCODE_TO_BCKCODE(SFI,CODE)\
        (GtCodetype) ((CODE) >> (SFI)->spmopt_kmerscancodeshift2bckcode)

#define GT_SCANCODE_TO_PREFIXCODE(SFI,CODE)\
        (GtCodetype) ((CODE) >> (SFI)->spmopt_kmerscancodeshift2prefixcode)

#if defined (_LP64) || defined (_WIN64)
#define GT_SCANCODE_TO_SUFFIXCODE(SFI,CODE)\
        (GtCodetype) ((CODE) & (SFI)->spmopt_kmerscancodesuffixmask)

static bool gt_checksuffixprefixbuckets(const Sfxiterator *sfi,
                                        GtCodetype scancode)
{
  GtCodetype prefixcode = GT_SCANCODE_TO_PREFIXCODE(sfi,scancode);
  GtCodetype suffixcode = GT_SCANCODE_TO_SUFFIXCODE(sfi,scancode);

  gt_assert(prefixcode < sfi->spmopt_numofallprefixcodes);
  gt_assert(suffixcode < sfi->spmopt_numofallsuffixcodes);
  return (GT_ISIBITSET(sfi->markprefixbuckets,prefixcode) &&
          GT_ISIBITSET(sfi->marksuffixbuckets,suffixcode)) ? true : false;
}
#else
static bool gt_checksuffixprefixbuckets(const Sfxiterator *sfi,
                                        GtCodetype scancode)
{
  GtCodetype prefixcode = GT_SCANCODE_TO_PREFIXCODE(sfi,scancode);

  gt_assert(prefixcode < sfi->spmopt_numofallprefixcodes);
  return GT_ISIBITSET(sfi->markprefixbuckets,prefixcode) ? true : false;
}
#endif

#define GT_INSERTKMERWITHOUTSPECIAL1(SFI,FIRSTINRANGE,POSITION,SEQNUM,RELPOS,\
                                     SCANCODE)\
        if ((SFI)->markprefixbuckets == NULL)\
        {\
          if ((SCANCODE) >= (SFI)->currentmincode &&\
              (SCANCODE) <= (SFI)->currentmaxcode)\
          {\
            GtUword stidx;\
            stidx = gt_bcktab_leftborder_insertionindex((SFI)->leftborder,\
                                                        SCANCODE);\
            /* from right to left */\
            GT_SUFFIXSORTSPACE_EXPORT_SET((SFI)->suffixsortspace,\
                                          (SFI)->exportptr,stidx,POSITION);\
          }\
        } else\
        {\
          GtCodetype bcktabcode = GT_SCANCODE_TO_BCKCODE((SFI),SCANCODE);\
          if (bcktabcode >= (SFI)->currentmincode &&\
              bcktabcode <= (SFI)->currentmaxcode &&\
              (FIRSTINRANGE || gt_checksuffixprefixbuckets(SFI,SCANCODE)))\
          {\
            GtUword stidx;\
            stidx = gt_bcktab_leftborder_insertionindex((SFI)->leftborder,\
                                                        bcktabcode);\
            /* from right to left */\
            GT_SUFFIXSORTSPACE_EXPORT_SET((SFI)->suffixsortspace,\
                                          (SFI)->exportptr,stidx,POSITION);\
          }\
        }

static void gt_insertkmerwithoutspecial(void *processinfo,
                                        GtUword position,
                                        const GtKmercode *kmercode)
{
  if (!kmercode->definedspecialposition)
  {
    GT_INSERTKMERWITHOUTSPECIAL1((Sfxiterator *) processinfo, false,
                                 position, 0, 0, kmercode->code);
  }
}

static void gt_reversespecialcodes(Codeatposition *spaceCodeatposition,
                                   GtUword nextfreeCodeatposition)
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
  GtUword insertindex, j, stidx;

  for (prefixindex=1U; prefixindex < sfi->prefixlength; prefixindex++)
  {
    for (j=0, insertindex = 0; j < sfi->nextfreeCodeatposition; j++)
    {
      if (prefixindex <= sfi->spaceCodeatposition[j].maxprefixindex)
      {
        code = gt_bcktab_codedownscale(sfi->bcktab,
                                       (GtCodetype)
                                       sfi->spaceCodeatposition[j].code,
                                       prefixindex,
                                       sfi->spaceCodeatposition[j].
                                            maxprefixindex);
        if (code >= sfi->currentmincode && code <= sfi->currentmaxcode)
        {
          gt_bcktab_updatespecials(sfi->bcktab,code,prefixindex);
          stidx = gt_bcktab_leftborder_insertionindex(sfi->leftborder,code);
          /* from right to left */
          gt_suffixsortspace_set(sfi->suffixsortspace,0,stidx,
                                 sfi->spaceCodeatposition[j].position
                                 - prefixindex);
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
  GtUword stidx;
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
            gt_bcktab_updatespecials(sfi->bcktab,code,prefixindex);
            gt_assert(code > 0);
            stidx = gt_bcktab_leftborder_insertionindex(sfi->leftborder,code);
            /* from right to left */
            gt_suffixsortspace_set(sfi->suffixsortspace,0,stidx,
                                   specialcontext.position - prefixindex);
          }
        }
      }
    }
    gt_Enumcodeatposition_delete(ecp);
    ecp = NULL;
  }
}

#ifdef GT_THREADS_ENABLED
#ifdef GT_THREADS_PARTITION
static void gt_suftabparts_multi_delete(GtSuftabparts **suftabparts_tab,
                                        unsigned int parts)
{
  if (suftabparts_tab != NULL)
  {
    unsigned int part;

    for (part = 0; part < parts; part++)
    {
      gt_suftabparts_delete(suftabparts_tab[part]);
    }
    gt_free(suftabparts_tab);
  }
}
#endif
#endif

int gt_Sfxiterator_delete(Sfxiterator *sfi,GtError *err)
{
  bool haserr = false;

  if (sfi == NULL)
  {
    return 0;
  }
#ifdef SKDEBUG
  if (sfi->bcktab != NULL)
  {
    gt_bcktab_checkcountspecialcodes(sfi->bcktab);
  }
#endif
  if (sfi->sri != NULL)
  {
    gt_specialrangeiterator_delete(sfi->sri);
  }
  gt_free(sfi->spaceCodeatposition);
  sfi->spaceCodeatposition = NULL;
  gt_suffixsortspace_delete(sfi->suffixsortspace,
                            sfi->sfxstrategy.spmopt_minlength == 0
                              ? true : false);
  if (sfi->suftabparts != NULL &&
      gt_suftabparts_numofparts(sfi->suftabparts) > 1U &&
      sfi->outfpbcktab != NULL)
  {
    gt_suftabparts_showallrecords(sfi->suftabparts,true);
    if (gt_bcktab_remap_all(sfi->bcktab,err) != 0)
    {
      haserr = true;
    } else
    {
      int ret = gt_bcktab_flush_to_file(sfi->outfpbcktab,sfi->bcktab,err);
      gt_fa_fclose(sfi->outfpbcktab);
      if (ret != 0)
      {
        haserr = true;
      }
    }
  }
  gt_bcktab_delete(sfi->bcktab);
#ifdef GT_THREADS_ENABLED
#ifdef GT_THREADS_PARTITION
  gt_suftabparts_multi_delete(sfi->partitions_for_threads,
                              gt_suftabparts_numofparts(sfi->suftabparts));
#endif
#endif
  gt_suftabparts_delete(sfi->suftabparts);
  gt_Outlcpinfo_delete(sfi->outlcpinfoforsample);
  if (sfi->mappedmarkprefixbuckets == NULL)
  {
    gt_free(sfi->markprefixbuckets);
  }
  gt_Sfxmappedrange_delete(sfi->mappedmarkprefixbuckets);
  sfi->mappedmarkprefixbuckets = NULL;
  gt_free(sfi->marksuffixbuckets);
  gt_differencecover_delete(sfi->dcov);
  gt_SSSPbuf_delete(sfi->sssp_buf);
  gt_free(sfi);
  return haserr ? -1 : 0;
}

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
    GtUword position = 0;

    while ((kmercodeptr = gt_kmercodeiterator_encseq_next(kmercodeiterator))
                          != NULL)
    {
      updatekmercount(sfi,position++,kmercodeptr);
    }
  }
  gt_kmercodeiterator_delete(kmercodeiterator);
}

void getencseqkmersinsertkmerwithoutspecial(const GtEncseq *encseq,
                                            GtReadmode readmode,
                                            unsigned int kmersize,
                                            Sfxiterator *sfi)
{
  GtKmercodeiterator *kmercodeiterator;
  const GtKmercode *kmercodeptr;

  kmercodeiterator = gt_kmercodeiterator_encseq_new(encseq,readmode,kmersize,0);
  if (!gt_kmercodeiterator_inputexhausted(kmercodeiterator))
  {
    GtUword position = 0;

    while ((kmercodeptr = gt_kmercodeiterator_encseq_next(kmercodeiterator))
                          != NULL)
    {
      gt_insertkmerwithoutspecial(sfi,position++,kmercodeptr);
    }
    gt_kmercodeiterator_delete(kmercodeiterator);
  }
}

#undef DEBUGSIZEESTIMATION
#ifdef DEBUGSIZEESTIMATION
static void verifyestimatedspace(size_t estimatedspace)
{
  GtUword usedspace_ma_fa = gt_ma_get_space_current() +
                                  gt_fa_get_space_current();
  if (usedspace_ma_fa > 0)
  {
    double relativedifference;

    if (usedspace_ma_fa >= (GtUword) estimatedspace)
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
      fprintf(stderr, "%s: relativedifference %.4f too large: "
                      "estimatedspace=%.4f, usedspace_ma_fa=%.4f\n",
                      __func__,relativedifference,
                      GT_MEGABYTES(estimatedspace),
                      GT_MEGABYTES(usedspace_ma_fa));
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
}
#endif

GtCodetype gt_kmercode_at_firstpos(const GtTwobitencoding *twobitencoding,
                                   unsigned int kmersize)
{
  const GtCodetype maskright = GT_MASKRIGHT(kmersize);
  return (GtCodetype) (twobitencoding[0] >>
                       GT_MULT2(GT_UNITSIN2BITENC - kmersize)) & maskright;
}

#ifdef SKDEBUG

static void checkallreversebitpairs(void)
{
  unsigned int kmersize, code, coderev, coderevrev, maxcode;

  for (kmersize = 2U; kmersize <= 14U; kmersize++)
  {
    maxcode = (1U << 2 * kmersize)-1;
    printf("kmsize=%u,maxcode=%u\n",kmersize,maxcode);
    for (code = 0; code <= maxcode; code++)
    {
      coderev = gt_kmercode_reverse(code,kmersize);
      coderevrev = gt_kmercode_reverse(coderev,kmersize);
      gt_assert(coderevrev != code);
    }
  }
}
#endif

#define GT_UPDATEKMER(KMER,CC,MASKRIGHT)\
        KMER <<= 2;\
        KMER |= CC;\
        KMER &= MASKRIGHT

#define GT_ADJUSTREVERSEPOS(RB,POS) ((RB) - (POS))

static GtCodetype getencseqkmers_nospecialtwobitencoding(
                                    const GtTwobitencoding *twobitencoding,
                                    GtUword totallength,
                                    GtCodetype maskright,
                                    GtReadmode readmode,
                                    unsigned int kmersize,
                                    unsigned int upperkmersize,
                                    void(*processkmercode)(void *,
                                                           bool,
                                                           GtUword,
                                                           GtCodetype),
                                    void *processkmercodeinfo,
                                    bool onlyfirst,
                                    GtUword startpos,
                                    GtUword endpos)
{
  GtUword pos, unitindex, rightbound = totallength - kmersize;
  unsigned int shiftright;
  GtCodetype code;
  GtUchar cc;
  GtTwobitencoding currentencoding;

  gt_assert(kmersize > 1U);
  if (GT_ISDIRREVERSE(readmode))
  {
    GtUword startpos2;

    gt_assert(endpos >= (GtUword) upperkmersize);
    pos = endpos - (GtUword) kmersize;
    unitindex = (pos > 0) ? GT_DIVBYUNITSIN2BITENC(pos-1) : 0;
    code = gt_kmercode_reverse(gt_kmercode_at_position(twobitencoding,pos,
                                                       kmersize),
                               kmersize);
    if (processkmercode != NULL)
    {
      processkmercode(processkmercodeinfo,
                      true,
                      GT_ADJUSTREVERSEPOS(rightbound,pos),
                      GT_ISDIRCOMPLEMENT(readmode)
                        ? gt_kmercode_complement(code,maskright)
                        : code);
    }
    if (onlyfirst)
    {
      return code;
    }
    currentencoding = twobitencoding[unitindex];
    startpos2 = startpos + (upperkmersize - kmersize);
    shiftright = (unsigned int)
                 GT_MULT2(GT_UNITSIN2BITENC - 1 -
                          GT_MODBYUNITSIN2BITENC(pos-1));
    while (pos > startpos2)
    {
      pos--;
      cc = (GtUchar) (currentencoding >> shiftright) & 3;
      GT_UPDATEKMER(code,cc,maskright);
      if (processkmercode != NULL)
      {
        processkmercode(processkmercodeinfo,false,
                        GT_ADJUSTREVERSEPOS(rightbound,pos),
                         (readmode == GT_READMODE_REVCOMPL)
                          ? gt_kmercode_complement(code,maskright)
                          : code);
      }
      if (shiftright < (unsigned int) (GT_INTWORDSIZE-2))
      {
        shiftright += 2;
      } else
      {
        gt_assert(unitindex > 0 || pos == startpos2);
        if (unitindex > 0)
        {
          currentencoding = twobitencoding[--unitindex];
        }
        shiftright = 0;
      }
    }
  } else
  {
    GtUword maxunitindex = gt_unitsoftwobitencoding(totallength) - 1;

    pos = startpos;
    unitindex = GT_DIVBYUNITSIN2BITENC(startpos+kmersize);
    code = gt_kmercode_at_position(twobitencoding,pos,kmersize);
    if (processkmercode != NULL)
    {
      processkmercode(processkmercodeinfo,true,pos,
                      GT_ISDIRCOMPLEMENT(readmode)
                        ? gt_kmercode_complement(code,maskright)
                        : code);
    }
    if (onlyfirst)
    {
      return code;
    }
    currentencoding = twobitencoding[unitindex];
    shiftright = (unsigned int)
                 GT_MULT2(GT_UNITSIN2BITENC - 1 -
                          GT_MODBYUNITSIN2BITENC(startpos+kmersize));
    while (pos < endpos - (GtUword) upperkmersize)
    {
      pos++;
      cc = (GtUchar) (currentencoding >> shiftright) & 3;
      GT_UPDATEKMER(code,cc,maskright);
      if (processkmercode != NULL)
      {
        processkmercode(processkmercodeinfo,false,pos,
                          (readmode == GT_READMODE_COMPL)
                            ? gt_kmercode_complement(code,maskright)
                            : code);
      }
      if (shiftright > 0)
      {
        shiftright -= 2;
      } else
      {
        gt_assert(unitindex < maxunitindex-1 ||
                  pos == endpos - (GtUword) upperkmersize);
        if (unitindex < maxunitindex-1)
        {
          currentencoding = twobitencoding[++unitindex];
        }
        shiftright = (unsigned int) (GT_INTWORDSIZE-2);
      }
    }
  }
  return code;
}

static
void getencseqkmers_rangetwobitencoding(const GtEncseq *encseq,
                                        GtReadmode readmode,
                                        unsigned int kmersize,
                                        unsigned int upperkmersize,
                                        bool onlyfirst,
                                        void(*processkmercode)(void *,
                                                               bool,
                                                               GtUword,
                                                               GtCodetype),
                                        void *processkmercodeinfo,
                                        void(*processkmerspecial)(void *,
                                                                  unsigned int,
                                                                  unsigned int,
                                                                  GtUword),
                                        void *processkmerspecialinfo,
                                        GtUword startpos,
                                        GtUword endpos)
{
  GtCodetype lastcode, newcode, realtotallength;
  const GtTwobitencoding *twobitencoding
          = gt_encseq_twobitencoding_export(encseq);
  const GtCodetype maskright = GT_MASKRIGHT(kmersize);
  bool mirrored = gt_encseq_is_mirrored(encseq);
  const GtUword totallength = gt_encseq_total_length(encseq);
  realtotallength = totallength;
  if (mirrored) {
    gt_assert((totallength & 1) == 1UL);
    realtotallength = ((totallength - 1) / 2);
  }

  if (mirrored && startpos >= realtotallength) {
    gt_readmode_invert(readmode);
    startpos = GT_REVERSEPOS(realtotallength, startpos - realtotallength - 2);
    if (endpos == totallength)
      endpos = 0;
    else
      endpos = GT_REVERSEPOS(realtotallength, endpos - realtotallength - 2);
    if (startpos > endpos) {
      GtUword tmp = startpos;
      startpos = endpos;
      endpos = tmp;
    }
    gt_assert(startpos <= endpos);
    gt_assert(endpos <= realtotallength);
  }
  if (endpos - startpos >= (GtUword) upperkmersize)
  {
    gt_assert(endpos > 0);
    lastcode = getencseqkmers_nospecialtwobitencoding(twobitencoding,
                                                      totallength,
                                                      maskright,
                                                      readmode,
                                                      kmersize,
                                                      upperkmersize,
                                                      processkmercode,
                                                      processkmercodeinfo,
                                                      onlyfirst,
                                                      startpos,
                                                      endpos);
    if (processkmerspecial != NULL)
    {
      if (GT_ISDIRCOMPLEMENT(readmode))
      {
        lastcode = gt_kmercode_complement(lastcode,maskright);
      }
      newcode = ((lastcode << 2) | 3UL) & maskright;
      processkmerspecial(processkmerspecialinfo,
                         kmersize-1,
                         (unsigned int) newcode,
                         GT_ISDIRREVERSE(readmode) ? (totallength - startpos)
                                                   : endpos);
    }
  } else
  {
    if (processkmerspecial != NULL && startpos < endpos)
    {
      unsigned int fillpos;

      gt_assert((GtUword) kmersize > endpos - startpos);
      fillpos = (unsigned int) (kmersize - (endpos - startpos));
      lastcode = gt_kmercode_at_position(twobitencoding,startpos,
                                         (unsigned int) (endpos - startpos));
      if (GT_ISDIRREVERSE(readmode) && (unsigned int) (endpos - startpos) > 1U)
      {
        lastcode = gt_kmercode_reverse(lastcode,
                                       (unsigned int) (endpos-startpos));
      }
      if (GT_ISDIRCOMPLEMENT(readmode))
      {
        lastcode = gt_kmercode_complement(lastcode,maskright);
      }
      newcode
        = ((lastcode << GT_MULT2(fillpos)) | ((1UL << GT_MULT2(fillpos)) - 1))
           & maskright;
      processkmerspecial(processkmerspecialinfo,
                         (unsigned int) (endpos - startpos),
                         (unsigned int) newcode,
                         GT_ISDIRREVERSE(readmode) ? (totallength - startpos)
                                                   : endpos);
    }
  }
}

void getencseqkmers_twobitencoding_slice(const GtEncseq *encseq,
                                         GtReadmode readmode,
                                         unsigned int kmersize,
                                         unsigned int upperkmersize,
                                         bool onlyfirst,
                                         void(*processkmercode)(void *,
                                                                bool,
                                                                GtUword,
                                                                GtCodetype),
                                         void *processkmercodeinfo,
                                         void(*processkmerspecial)(void *,
                                                                   unsigned int,
                                                                   unsigned int,
                                                                   GtUword),
                                         void *processkmerspecialinfo,
                                         GtUword slice_startpos,
                                         GtUword slice_endpos)
{
  GtUword laststart = slice_startpos,
          lastend = slice_endpos;
  gt_assert(laststart <= lastend);

  if (gt_encseq_has_specialranges(encseq))
  {
    GtSpecialrangeiterator *sri;
    GtRange range;

    if (GT_ISDIRREVERSE(readmode))
    {
      sri = gt_specialrangeiterator_new(encseq,false);
      while (gt_specialrangeiterator_next(sri,&range))
      {
        if (range.end <= lastend) {
          getencseqkmers_rangetwobitencoding(encseq,
                                             readmode,
                                             kmersize,
                                             upperkmersize,
                                             onlyfirst,
                                             processkmercode,
                                             processkmercodeinfo,
                                             processkmerspecial,
                                             processkmerspecialinfo,
                                             range.end,
                                             lastend);
          lastend = range.start;
        }
        if (laststart > range.start) {
          break;
        }
      }
    } else
    {
      sri = gt_specialrangeiterator_new(encseq,true);
      while (gt_specialrangeiterator_next(sri,&range))
      {
        if (range.start >= laststart) {
          getencseqkmers_rangetwobitencoding(encseq,
                                             readmode,
                                             kmersize,
                                             upperkmersize,
                                             onlyfirst,
                                             processkmercode,
                                             processkmercodeinfo,
                                             processkmerspecial,
                                             processkmerspecialinfo,
                                             laststart,
                                             range.start);
          laststart = range.end;
        }
        if (lastend < range.end) {
          break;
        }
      }
    }
    gt_assert(gt_encseq_total_length(encseq) >= laststart);
    gt_specialrangeiterator_delete(sri);
  }

  if (laststart <= slice_endpos && slice_startpos <= lastend) {
    GtUword startpos = GT_ISDIRREVERSE(readmode) ? slice_startpos : laststart;
    GtUword endpos = GT_ISDIRREVERSE(readmode) ? lastend : slice_endpos;

    getencseqkmers_rangetwobitencoding(encseq,
                                       readmode,
                                       kmersize,
                                       upperkmersize,
                                       onlyfirst,
                                       processkmercode,
                                       processkmercodeinfo,
                                       processkmerspecial,
                                       processkmerspecialinfo,
                                       startpos,
                                       endpos);
  }
}

void getencseqkmers_twobitencoding(const GtEncseq *encseq,
                                   GtReadmode readmode,
                                   unsigned int kmersize,
                                   unsigned int upperkmersize,
                                   bool onlyfirst,
                                   void(*processkmercode)(void *,
                                                          bool,
                                                          GtUword,
                                                          GtCodetype),
                                   void *processkmercodeinfo,
                                   void(*processkmerspecial)(void *,
                                                             unsigned int,
                                                             unsigned int,
                                                             GtUword),
                                   void *processkmerspecialinfo)
{
  getencseqkmers_twobitencoding_slice(encseq,
                                      readmode,
                                      kmersize,
                                      upperkmersize,
                                      onlyfirst,
                                      processkmercode,
                                      processkmercodeinfo,
                                      processkmerspecial,
                                      processkmerspecialinfo,
                                      0,
                                      gt_encseq_total_length(encseq));
}

static void gt_updateleftborderforkmer(Sfxiterator *sfi,
                                       GT_UNUSED bool firstinrange,
                                       GT_UNUSED GtUword position,
                                       GT_UNUSED GtUword seqnum,
                                       GT_UNUSED GtUword relpos,
                                       GtCodetype code)
{
  gt_assert(sfi->sfxstrategy.spmopt_minlength == 0);
  gt_bcktab_leftborder_addcode(sfi->leftborder,code);
}

static void gt_updateleftborderforspecialkmer(Sfxiterator *sfi,
                                              unsigned int maxprefixindex,
                                              GtUword position,
                                              unsigned int code)
{
  unsigned int idx;

  gt_assert(sfi->sfxstrategy.spmopt_minlength == 0);
  if (sfi->sfxstrategy.storespecialcodes)
  {
    Codeatposition *spcaptr;
    spcaptr = sfi->spaceCodeatposition + sfi->nextfreeCodeatposition++;
    spcaptr->maxprefixindex = maxprefixindex;
    spcaptr->code = code;
    spcaptr->position = position;
  }
  for (idx=maxprefixindex; idx>=1U; idx--)
  {
    gt_bcktab_leftborder_addcode(sfi->leftborder,(GtCodetype) code);
    code = ((code << 2) | 3U) & sfi->kmerfastmaskright;
  }
}

#define GT_SPMOPT_UPDATELEFTBORDERFORKMER(SFI,FIRSTINRANGE,POSITION,SEQNUM,\
                                          RELPOS,SCANCODE)\
        gt_assert((SFI)->sfxstrategy.spmopt_minlength > 0);\
        if (FIRSTINRANGE || gt_checksuffixprefixbuckets(SFI,SCANCODE))\
        {\
          gt_bcktab_leftborder_addcode((SFI)->leftborder,\
                                       GT_SCANCODE_TO_BCKCODE(SFI,SCANCODE));\
        }

typedef struct
{
  const GtTwobitencoding *twobitencoding;
  GtUword totallength, maxunitindex, realtotallength, rightbound,
                numofsequences;
  GtCodetype maskright;
  unsigned int kmersize, upperkmersize;
  bool mirrored;
  const GtEncseq *encseq; /* XXX remove later */
} GtSfxmapped4constinfo;

/* This is for checking only
#define GT_ENCSEQ_RELPOS_SEQNUM_CHECK(POS)\
        gt_encseq_relpos_seqnum_check(__FILE__,__LINE__,\
                                      mapped4info->encseq,relpos,\
                                      specialfreeunit,POS)
*/
#define GT_ENCSEQ_RELPOS_SEQNUM_CHECK(POS) /* Nothing */

#define PROCESSKMERPREFIX(FUN) updateleftborder_##FUN
#define PROCESSKMERTYPE        Sfxiterator
#define PROCESSKMERSPECIALTYPE GT_UNUSED Sfxiterator
#define PROCESSKMERCODE        gt_updateleftborderforkmer
#define PROCESSKMERCODESPECIAL gt_updateleftborderforspecialkmer

#include "sfx-mapped4.gen"

#undef PROCESSKMERPREFIX
#undef PROCESSKMERTYPE
#undef PROCESSKMERSPECIALTYPE
#undef PROCESSKMERCODE
#undef PROCESSKMERCODESPECIAL

/* start with next inling */

#define PROCESSKMERPREFIX(FUN) spmopt_updateleftborder_##FUN
#define PROCESSKMERTYPE        Sfxiterator
#define PROCESSKMERSPECIALTYPE GT_UNUSED Sfxiterator
#define PROCESSKMERCODE        GT_SPMOPT_UPDATELEFTBORDERFORKMER
#define GT_IGNORERIGHTBOUND

#include "sfx-mapped4.gen"
#undef GT_IGNORERIGHTBOUND

/* start with next inling */

#undef PROCESSKMERPREFIX
#undef PROCESSKMERTYPE
#undef PROCESSKMERSPECIALTYPE
#undef PROCESSKMERCODE

#define PROCESSKMERPREFIX(FUN)          insertsuffix_##FUN
#define PROCESSKMERTYPE                 Sfxiterator
#define PROCESSKMERSPECIALTYPE          GT_UNUSED Sfxiterator
#define PROCESSKMERCODE                 GT_INSERTKMERWITHOUTSPECIAL1

#include "sfx-mapped4.gen"

#undef PROCESSKMERPREFIX
#undef PROCESSKMERTYPE
#undef PROCESSKMERSPECIALTYPE
#undef PROCESSKMERCODE

/*
#define SHOWCURRENTSPACE\
        printf("spacepeak at line %d: %.2f\n",__LINE__,\
          GT_MEGABYTES(gt_ma_get_space_current() + gt_fa_get_space_current()))
*/

#define SHOWCURRENTSPACE /* Nothing */

static void gt_sfimarkprefixsuffixbuckets(void *processinfo,
                                          GT_UNUSED bool firstinrange,
                                          GT_UNUSED GtUword pos,
                                          GtCodetype scancode)
{
  Sfxiterator *sfi = (Sfxiterator *) processinfo;
  GtCodetype checkcode = GT_SCANCODE_TO_PREFIXCODE(sfi,scancode);

  gt_assert(firstinrange);
  if (!GT_ISIBITSET(sfi->markprefixbuckets,checkcode))
  {
    GT_SETIBIT(sfi->markprefixbuckets,checkcode);
  }
#if defined (_LP64) || defined (_WIN64)
  checkcode = GT_SCANCODE_TO_SUFFIXCODE(sfi,scancode);
  if (!GT_ISIBITSET(sfi->marksuffixbuckets,checkcode))
  {
    GT_SETIBIT(sfi->marksuffixbuckets,checkcode);
  }
#endif
}

static size_t gt_sizeforbittable(unsigned int numofchars,
                                 unsigned int prefixlength)
{
  GtUword numofcodes;

  numofcodes = gt_power_for_small_exponents(numofchars,prefixlength);
  return sizeof (GtBitsequence) * GT_NUMOFINTSFORBITS(numofcodes);
}

static void gt_determineaddionalsuffixprefixchars(
                                   unsigned int *additionalprefixchars,
                                   unsigned int *additionalsuffixchars,
                                   unsigned int numofchars,
                                   unsigned int prefixlength,
                                   size_t estimatedspace,
                                   GtUword maximumspace)
{
  unsigned int prefixchars;
  size_t sizeofprefixmarks;

  for (prefixchars = 1U;
       prefixlength + prefixchars <= (unsigned int) GT_UNITSIN2BITENC;
       prefixchars++)
  {
    sizeofprefixmarks = gt_sizeforbittable(numofchars,prefixlength+prefixchars);
    if (estimatedspace + sizeofprefixmarks > (size_t) maximumspace)
    {
      prefixchars--;
      break;
    }
  }
  *additionalprefixchars = prefixchars;
#if defined (_LP64) || defined (_WIN64)
  {
    unsigned int suffixchars;
    size_t sizeofsuffixmarks;
    sizeofprefixmarks = gt_sizeforbittable(numofchars,prefixlength+prefixchars);
    for (suffixchars = 1U;
         prefixlength+prefixchars+prefixlength+suffixchars <=
           (unsigned) GT_UNITSIN2BITENC;
         suffixchars++)
    {
      sizeofsuffixmarks = gt_sizeforbittable(numofchars,
                                             prefixlength+suffixchars);
      if (estimatedspace + sizeofprefixmarks + sizeofsuffixmarks
          > (size_t) maximumspace)
      {
        suffixchars--;
        break;
      }
    }
    if (prefixchars <= suffixchars)
    {
      suffixchars = prefixchars - 1;
    }
    *additionalsuffixchars = suffixchars;
  }
#else
  *additionalsuffixchars = 0;
#endif
}

static GtUword gt_bcktab_code_to_prefix_index(GtUword code,
                                                    const void *data)
{
  unsigned int additionalprefixchars, *ptr = (unsigned int *) data;

  additionalprefixchars = *ptr;
  if (GT_MULT2(additionalprefixchars) > (unsigned int) GT_LOGWORDSIZE)
  {
    return (GtUword) (code << (GT_MULT2(additionalprefixchars) -
                                 GT_LOGWORDSIZE));
  }
  return (GtUword) (code >> (GT_LOGWORDSIZE -
                                GT_MULT2(additionalprefixchars)));
}

static void gt_bcktab_code_to_minmax_prefix_index(GtUword *mincode,
                                                  GtUword *maxcode,
                                                  const void *data)
{
  *mincode = gt_bcktab_code_to_prefix_index(*mincode,data);
  *maxcode = gt_bcktab_code_to_prefix_index(*maxcode,data);
}

#ifdef GT_THREADS_ENABLED
#define GT_SFX_THREADS_JOBS gt_jobs
#ifdef GT_THREADS_PARTITION
static GtSuftabparts **gt_partitions_for_threads_new(
                               const GtSuftabparts *suftabparts,
                               const GtBcktab *bcktab,
                               GtSfxmappedrangelist *sfxmrlist,
                               GtUword specialcharacters,
                               GtLogger *logger)
{
  if (GT_SFX_THREADS_JOBS > 1U)
  {
    const unsigned int parts = gt_suftabparts_numofparts(suftabparts);
    unsigned int part;
    GtSuftabparts **partitions_for_threads;

    gt_assert(parts > 0);
    partitions_for_threads = gt_malloc(sizeof *partitions_for_threads * parts);
    /*gt_bcktab_leftborder_show(bcktab);*/
    gt_suftabparts_showallrecords(suftabparts,true);
    for (part = 0; part < parts; part++)
    {
      unsigned int parts_sub;
      GtCodetype mincode = gt_suftabparts_minindex(part,suftabparts);
      GtCodetype maxcode = gt_suftabparts_maxindex(part,suftabparts);
      GtUword widthofpart = gt_suftabparts_widthofpart(part,suftabparts);

      partitions_for_threads[part]
        = gt_suftabparts_new(GT_SFX_THREADS_JOBS,
                             bcktab,
                             mincode,
                             maxcode,
                             NULL,
                             sfxmrlist,
                             widthofpart,
                             specialcharacters + 1,
                             logger);
      gt_suftabparts_showallrecords(partitions_for_threads[part],true);
      parts_sub = gt_suftabparts_numofparts(partitions_for_threads[part]);
      gt_assert(parts_sub > 0);
      if (mincode != gt_suftabparts_minindex(0,partitions_for_threads[part]))
      {
        fprintf(stderr,"part %u: mincode=" GT_WU "!=" GT_WU "=minindex[0]\n",
                        part,
                        mincode,
                        gt_suftabparts_minindex(0,partitions_for_threads[part]))
                        ;
        exit(EXIT_FAILURE);
      }
      if (maxcode != gt_suftabparts_maxindex(parts_sub - 1,
                                             partitions_for_threads[part]))
      {
        fprintf(stderr,"part %u: maxcode=" GT_WU " != " GT_WU "maxindex[%u]\n",
                        part,
                        maxcode,
                        gt_suftabparts_maxindex(parts_sub-1,
                                                partitions_for_threads[part]),
                        parts_sub-1);
        exit(EXIT_FAILURE);
      }
    }
    return partitions_for_threads;
  } else
  {
    return NULL;
  }
}
#endif
#else
#define GT_SFX_THREADS_JOBS 1U
#endif

Sfxiterator *gt_Sfxiterator_new_withadditionalvalues(
                                const GtEncseq *encseq,
                                GtReadmode readmode,
                                unsigned int prefixlength,
                                unsigned int numofparts,
                                GtUword maximumspace,
                                void *voidoutlcpinfo,
                                FILE *outfpbcktab,
                                const Sfxstrategy *sfxstrategy,
                                GtTimer *sfxprogress,
                                bool withprogressbar,
                                GtLogger *logger,
                                GtError *err)
{
  Sfxiterator *sfi = NULL;
  GtUword realspecialranges, specialcharacters, numofsuffixestosort = 0;
  bool haserr = false;
  GtSfxmappedrangelist *sfxmrlist = gt_Sfxmappedrangelist_new();
#if defined (_LP64) || defined (_WIN64)
  size_t estimatedspace = (size_t) 13131;
#else
  size_t estimatedspace = (size_t) 7968;
#endif

  gt_error_check(err);
  SHOWCURRENTSPACE;
  gt_assert(encseq != NULL);
  estimatedspace += (size_t) gt_encseq_sizeofrep(encseq) +
                             gt_encseq_sizeofstructure();
  realspecialranges = gt_encseq_realspecialranges(encseq);
  specialcharacters = gt_encseq_specialcharacters(encseq);
  gt_assert(prefixlength > 0);
  if (sfxstrategy != NULL)
  {
    if (sfxstrategy->storespecialcodes &&
        prefixlength > (unsigned int) GT_MAXPREFIXLENGTH)
    {
      gt_error_set(err,"argument for option -pl must be in the range [1,%u]",
                   GT_MAXPREFIXLENGTH);
      haserr = true;
    } else
    {
      if (sfxstrategy->spmopt_minlength > 0 &&
          prefixlength > sfxstrategy->spmopt_minlength)
      {
        gt_error_set(err,"argument for option -pl must not be larger "
                         "than argument to option -spmopt");
        haserr = true;
      }
    }
  }
  if (!haserr)
  {
    sfi = gt_malloc(sizeof (*sfi));
    estimatedspace += sizeof (*sfi);
    if (sfxstrategy != NULL && sfxstrategy->storespecialcodes &&
        sfxstrategy->spmopt_minlength == 0)
    {
      sfi->spaceCodeatposition
        = gt_malloc(sizeof (*sfi->spaceCodeatposition) * (realspecialranges+1));
      gt_logger_log(logger,"sizeof (spaceCodeatposition)="GT_WU" bytes",
                 (GtUword) (sizeof (*sfi->spaceCodeatposition) *
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
#ifdef GT_THREADS_ENABLED
#ifdef GT_THREADS_PARTITION
    sfi->partitions_for_threads = NULL;
#endif
#endif
    sfi->encseq = encseq;
    sfi->readmode = readmode;
    sfi->numofchars = gt_encseq_alphabetnumofchars(encseq);
    sfi->prefixlength = prefixlength;
    sfi->kmerfastmaskright = (1U << GT_MULT2(prefixlength))-1;
    sfi->mappedmarkprefixbuckets = NULL;
    sfi->markprefixbuckets = NULL;
    sfi->marksuffixbuckets = NULL;
    sfi->outfpbcktab = outfpbcktab;
    sfi->spmopt_kmerscansize = 0;
    sfi->spmopt_numofallprefixcodes = 0;
    sfi->spmopt_numofallsuffixcodes = 0;
    sfi->spmopt_kmerscancodeshift2bckcode = 0;
    sfi->spmopt_kmerscancodeshift2prefixcode = 0;
    sfi->spmopt_kmerscancodesuffixmask = 0;
    sfi->spmopt_additionalprefixchars = 3U;
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
    gt_logger_log(logger,"maxinsertionsort="GT_WU"",
                  sfi->sfxstrategy.maxinsertionsort);
    gt_logger_log(logger,"maxbltriesort="GT_WU"",
                  sfi->sfxstrategy.maxbltriesort);
    gt_logger_log(logger,"maxcountingsort="GT_WU"",
                  sfi->sfxstrategy.maxcountingsort);
    gt_logger_log(logger,"storespecialcodes=%s",
                  sfi->sfxstrategy.storespecialcodes ? "true" : "false");
    gt_logger_log(logger,"cmpcharbychar=%s",
                  sfi->sfxstrategy.cmpcharbychar ? "true" : "false");
    sfi->totallength = gt_encseq_total_length(encseq);
    gt_logger_log(logger,"totallength="GT_WU"",sfi->totallength);
    sfi->specialcharacters = specialcharacters;
    sfi->outlcpinfo = (GtOutlcpinfo *) voidoutlcpinfo;
    sfi->outlcpinfoforsample = NULL;
    sfi->sri = NULL;
    sfi->part = 0;
    sfi->exhausted = false;
    sfi->bucketiterstep = 0;
    sfi->logger = logger;
    sfi->sssp_buf = NULL;
    sfi->sfxprogress = sfxprogress;
    if (sfi->sfxstrategy.differencecover > 0 &&
        specialcharacters < sfi->totallength)
    {
      if (sfi->outlcpinfo != NULL)
      {
        sfi->outlcpinfoforsample
          = gt_Outlcpinfo_new(NULL,sfi->numofchars,0,false,false,NULL,NULL,err);
        if (sfi->outlcpinfoforsample == NULL)
        {
          haserr = true;
        }
      }
      if (!haserr)
      {
        /* the following function only has an effect for differencecover > 0 */
        sfi->dcov = gt_differencecover_prepare_sample(
                                     sfi->sfxstrategy.differencecover,
                                     sfi->encseq,
                                     sfi->readmode,
                                     sfi->sfxstrategy.samplewithprefixlengthnull
                                       ? 0 : sfi->prefixlength,
                                     &sfi->sfxstrategy,
                                     sfi->outlcpinfoforsample,
                                     sfi->logger,
                                     sfi->sfxprogress,
                                     err);
        if (sfi->dcov == NULL)
        {
          haserr = true;
        } else
        {
          if (gt_differencecover_is_empty(sfi->dcov))
          {
            gt_differencecover_delete(sfi->dcov);
            sfi->dcov = NULL;
          } else
          {
            size_t dcovspace = gt_differencecover_requiredspace(sfi->dcov);
            gt_logger_log(sfi->logger,"difference cover%srequires %.2f MB"
                           " (%.2f bytes/sample position)",
                          sfi->outlcpinfoforsample != NULL
                            ? " (including RMQ) "
                            : " ",
                          GT_MEGABYTES(dcovspace),
                          (double) dcovspace/
                                   gt_differencecover_samplesize(sfi->dcov));
            estimatedspace += dcovspace;
          }
        }
      }
    }
  }
  if (!haserr)
  {
    bool withspecialsuffixes;
    gt_assert(sfi != NULL);
    withspecialsuffixes = sfi->sfxstrategy.spmopt_minlength == 0 ? true : false;
    sfi->bcktab = gt_bcktab_new(sfi->numofchars,
                                prefixlength,
                                sfi->totallength+1,
                                sfi->sfxstrategy.storespecialcodes,
                                withspecialsuffixes,
                                sfi->logger,
                                err);
    if (sfi->bcktab == NULL)
    {
      sfi->leftborder = NULL;
      haserr = true;
    } else
    {
      uint64_t sizeofbcktab;

      sfi->leftborder = gt_bcktab_leftborder(sfi->bcktab);
      sizeofbcktab
        = gt_bcktab_sizeoftable(sfi->numofchars,prefixlength,
                                sfi->totallength+1,
                                sfi->sfxstrategy.spmopt_minlength == 0 ? true
                                                                       : false);
      estimatedspace += (size_t) sizeofbcktab +
                        gt_bcktab_sizeofworkspace(prefixlength);
    }
    SHOWCURRENTSPACE;
    if (prefixlength > 1U && gt_encseq_has_twobitencoding(sfi->encseq) &&
        sfi->sfxstrategy.spmopt_minlength > 0)
    {
      unsigned int suffixchars = 0, additionalsuffixchars = 2U;
      size_t sizeofprefixmarks, intsforbits;
#if defined (_LP64) || defined (_WIN64)
      size_t sizeofsuffixmarks;
#endif
      if (maximumspace > 0)
      {
        gt_determineaddionalsuffixprefixchars(
                               &sfi->spmopt_additionalprefixchars,
                               &additionalsuffixchars,
                               sfi->numofchars,
                               prefixlength,
                               estimatedspace,
                               maximumspace);
      }
#if defined (_LP64) || defined (_WIN64)
      if (sfi->prefixlength + sfi->spmopt_additionalprefixchars +
          sfi->prefixlength + additionalsuffixchars >
          (unsigned int) GT_UNITSIN2BITENC)
      {
        suffixchars = (unsigned int) (GT_UNITSIN2BITENC -
                      (sfi->prefixlength + sfi->spmopt_additionalprefixchars));
      } else
      {
        suffixchars = sfi->prefixlength + additionalsuffixchars;
      }
#endif
      sfi->spmopt_kmerscansize = sfi->prefixlength +
                                 sfi->spmopt_additionalprefixchars +
                                 suffixchars;
      gt_assert(sfi->spmopt_kmerscansize <= (unsigned int) GT_UNITSIN2BITENC);
      sfi->spmopt_kmerscancodeshift2bckcode
        = GT_MULT2(sfi->spmopt_additionalprefixchars + suffixchars);
      sfi->spmopt_kmerscancodeshift2prefixcode
        = GT_MULT2(suffixchars);
      sfi->spmopt_kmerscancodesuffixmask
        = (GtCodetype) ((1UL << GT_MULT2(suffixchars)) - 1);
      sfi->spmopt_numofallprefixcodes
        = gt_power_for_small_exponents(sfi->numofchars,
                                       sfi->prefixlength +
                                       sfi->spmopt_additionalprefixchars);
      GT_INITBITTAB(sfi->markprefixbuckets,
                    sfi->spmopt_numofallprefixcodes);
      intsforbits = GT_NUMOFINTSFORBITS(sfi->spmopt_numofallprefixcodes);
      sizeofprefixmarks = sizeof (*sfi->markprefixbuckets) * intsforbits;
      estimatedspace += sizeofprefixmarks;
      gt_logger_log(sfi->logger,"for all sequences, keep track of "
                                "%u-mers starting at position 0 using a "
                                "table of "GT_WU" bytes",
                                sfi->prefixlength +
                                sfi->spmopt_additionalprefixchars,
                                (GtUword) sizeofprefixmarks);
      sfi->mappedmarkprefixbuckets
        = gt_Sfxmappedrange_new("markprefixbuckets",
                                sfi->spmopt_numofallprefixcodes,
                                GtSfxGtBitsequence,
                                gt_bcktab_code_to_minmax_prefix_index,
                                &sfi->spmopt_additionalprefixchars);
      gt_Sfxmappedrangelist_add(sfxmrlist,sfi->mappedmarkprefixbuckets);
      sfi->spmopt_numofallsuffixcodes
        = gt_power_for_small_exponents(sfi->numofchars,suffixchars);
#if defined (_LP64) || defined (_WIN64)
      GT_INITBITTAB(sfi->marksuffixbuckets,
                    sfi->spmopt_numofallsuffixcodes);
      sizeofsuffixmarks
        = sizeof (*sfi->marksuffixbuckets) *
                  GT_NUMOFINTSFORBITS(sfi->spmopt_numofallsuffixcodes);
      estimatedspace += sizeofsuffixmarks;
      gt_logger_log(sfi->logger,"for all sequences, keep track of "
                                "%u-mers starting at position %u "
                                "using a table of "GT_WU" bytes",
                                suffixchars,
                                sfi->prefixlength +
                                sfi->spmopt_additionalprefixchars,
                                (GtUword) sizeofsuffixmarks);
#endif
      getencseqkmers_twobitencoding(encseq,
                                    sfi->readmode,
                                    sfi->spmopt_kmerscansize,
                                    sfi->spmopt_kmerscansize,
                                    true,
                                    gt_sfimarkprefixsuffixbuckets,
                                    sfi,
                                    NULL,
                                    NULL);
      if (maximumspace > 0)
      {
        gt_assert(estimatedspace <= (size_t) maximumspace);
      }
      /*printf("estimated space %.2f\n",GT_MEGABYTES(estimatedspace));*/
    }
  }
  SHOWCURRENTSPACE;
  if (!haserr)
  {
    GtUword largestbucketsize, saved_bucketswithoutwholeleaf;
    gt_assert(sfi != NULL);
    sfi->storespecials = true;
    if (sfxprogress != NULL)
    {
      gt_timer_show_progress(sfxprogress,"counting prefix distribution",stdout);
    }
    if (prefixlength == 1U)
    {
      unsigned int charidx;

      for (charidx=0; charidx<sfi->numofchars; charidx++)
      {
        unsigned int updateindex = GT_ISDIRCOMPLEMENT(readmode)
                                     ? GT_COMPLEMENTBASE(charidx)
                                     : charidx;
        gt_bcktab_leftborder_assign(sfi->leftborder,(GtCodetype) updateindex,
                                    gt_encseq_charcount(encseq,
                                                        (GtUchar) charidx));
      }
    } else
    {
      if (gt_encseq_has_twobitencoding(encseq) &&
          !sfi->sfxstrategy.kmerswithencseqreader)
      {
        if (sfi->sfxstrategy.spmopt_minlength == 0)
        {
          updateleftborder_getencseqkmers_twobitencoding(encseq,
                                                         readmode,
                                                         prefixlength,
                                                         prefixlength,
                                                         sfi,
                                                         sfi);
        } else
        {
          gt_assert(sfi->spmopt_kmerscansize > prefixlength);
          spmopt_updateleftborder_getencseqkmers_twobitencoding(
                                            encseq,
                                            readmode,
                                            sfi->spmopt_kmerscansize,
                                            sfi->sfxstrategy.spmopt_minlength,
                                            sfi,
                                            NULL);
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
      }
      if (sfi->sfxstrategy.storespecialcodes &&
          sfi->sfxstrategy.spmopt_minlength == 0)
      {
        size_t size_codeatposition = sizeof (*sfi->spaceCodeatposition) *
                                     sfi->nextfreeCodeatposition;
        gt_assert(realspecialranges+1 >= sfi->nextfreeCodeatposition);
        gt_logger_log(sfi->logger, "size for Codeatposition: %.2f MB=%.2f "
                                   "bytes/special suffix",
                 GT_MEGABYTES(size_codeatposition),
                 (double) size_codeatposition/specialcharacters);
        gt_reversespecialcodes(sfi->spaceCodeatposition,
                               sfi->nextfreeCodeatposition);
#ifdef SKDEBUG
        verifycodelistcomputation(encseq,
                                  readmode,
                                  realspecialranges,
                                  prefixlength,
                                  sfi->numofchars,
                                  sfi->nextfreeCodeatposition,
                                  sfi->spaceCodeatposition);
#endif
      }
    }
#ifdef SKDEBUG
    gt_bcktab_leftborder_show(sfi->bcktab);
#endif
    largestbucketsize
      = gt_bcktab_leftborderpartialsums(&saved_bucketswithoutwholeleaf,
                                        &numofsuffixestosort,
                                        sfi->bcktab);
    gt_logger_log(sfi->logger, "largest bucket size="GT_WU"",largestbucketsize);
    if (sfi->sfxstrategy.spmopt_minlength > 0)
    {
      gt_logger_log(sfi->logger, "relevant suffixes=%.2f%%",100.0 *
                                        (double) numofsuffixestosort/
                                        (sfi->totallength+1));
      gt_logger_log(sfi->logger,"saved_bucketswithoutwholeleaf="GT_WU"",
                                   saved_bucketswithoutwholeleaf);
      /*
      gt_assert(saved_bucketswithoutwholeleaf +
                numofsuffixestosort ==
                sfi->totallength - specialcharacters);
      */
    }
    if (sfi->outlcpinfo != NULL)
    {
      gt_Outlcpinfo_numsuffixes2output_set(
                                      sfi->outlcpinfo,
                                      sfi->sfxstrategy.spmopt_minlength == 0
                                        ? sfi->totallength + 1
                                        : numofsuffixestosort);
    }
    estimatedspace += sizeof (uint8_t) * largestbucketsize;
    SHOWCURRENTSPACE;
#ifdef DEBUGSIZEESTIMATION
    if (sfi->sfxstrategy.outsuftabonfile)
    {
      verifyestimatedspace(estimatedspace);
    }
#endif
    gt_bcktab_maprange_lb_cs(sfxmrlist,sfi->bcktab);
    if (maximumspace > 0)
    {
      int retval;

      gt_assert(numofparts == 1U);
      retval = gt_suftabparts_fit_memlimit(estimatedspace,
                                           maximumspace,
                                           sfi->bcktab,
                                           NULL,
                                           sfxmrlist,
                                           sfi->totallength,
                                           0, /* bitsforseqnumrelpos not
                                                 needed */
                                           specialcharacters,
                                           numofsuffixestosort,
                                           sfi->sfxstrategy.suftabuint,
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
/*
#define SHOWACTUALSPACE printf("line %d,realspace=heap=%.2f,map=%.2f\n",\
                              __LINE__,\
                              GT_MEGABYTES(gt_ma_get_space_current()),\
                              GT_MEGABYTES(gt_fa_get_space_current()))
*/
#define SHOWACTUALSPACE /* Nothing */
  SHOWACTUALSPACE;
  if (!haserr)
  {
    gt_assert(sfi != NULL);
    sfi->suftabparts = gt_suftabparts_new(numofparts,
                                         sfi->bcktab,
                                         (GtCodetype) 1,
                                         (GtCodetype) 0,
                                         NULL,
                                         sfxmrlist,
                                         numofsuffixestosort,
                                         specialcharacters + 1,
                                         logger);
    gt_assert(sfi->suftabparts != NULL);
#ifdef GT_THREADS_ENABLED
#ifdef GT_THREADS_PARTITION
    if (gt_suftabparts_numofparts(sfi->suftabparts) > 0)
    {
      sfi->partitions_for_threads
         = gt_partitions_for_threads_new(sfi->suftabparts,
                                         sfi->bcktab,
                                         sfxmrlist,
                                         specialcharacters,
                                         logger);
    }
#endif
#endif
    if (gt_suftabparts_numofparts(sfi->suftabparts) > 1U)
    {
      gt_bcktab_storetmp(sfi->bcktab);
    }
    SHOWACTUALSPACE;
    gt_assert(sfi != NULL && sfi->suftabparts != NULL);
    if (gt_suftabparts_numofparts(sfi->suftabparts) > 1U &&
        sfi->sfxstrategy.spmopt_minlength > 0)
    {
      gt_assert(sfi->markprefixbuckets != NULL);
      gt_assert(sfi->mappedmarkprefixbuckets != NULL);
      gt_Sfxmappedrange_storetmp_bitsequence(sfi->mappedmarkprefixbuckets,
                                             &sfi->markprefixbuckets,
                                             false);
      gt_assert(sfi->markprefixbuckets == NULL);
    } else
    {
      gt_Sfxmappedrange_delete(sfi->mappedmarkprefixbuckets);
      sfi->mappedmarkprefixbuckets = NULL;
    }
  }
  SHOWACTUALSPACE;
  if (!haserr)
  {
    gt_assert(sfi != NULL);
    if (gt_encseq_has_specialranges(sfi->encseq))
    {
      sfi->sri = gt_specialrangeiterator_new(sfi->encseq,
                                             GT_ISDIRREVERSE(sfi->readmode)
                                               ? false : true);
    } else
    {
      sfi->sri = NULL;
    }
    sfi->suffixsortspace
      = gt_suffixsortspace_new(gt_suftabparts_largest_width(sfi->suftabparts),
                               sfi->totallength,
                               sfi->sfxstrategy.suftabuint,
                               logger);
    gt_assert(sfi->suffixsortspace != NULL);
    sfi->sssp_buf
      = gt_SSSPbuf_new(gt_suftabparts_largest_width(sfi->suftabparts));
  }
  SHOWACTUALSPACE;
  gt_Sfxmappedrangelist_delete(sfxmrlist);
  if (haserr)
  {
    (void) gt_Sfxiterator_delete(sfi,NULL);
    return NULL;
  }
  return sfi;
}

Sfxiterator *gt_Sfxiterator_new(const GtEncseq *encseq,
                                GtReadmode readmode,
                                unsigned int prefixlength,
                                unsigned int numofparts,
                                GtUword maximumspace,
                                const Sfxstrategy *sfxstrategy,
                                GtTimer *sfxprogress,
                                bool withprogressbar,
                                GtLogger *logger,
                                GtError *err)
{
  return gt_Sfxiterator_new_withadditionalvalues(
                                encseq,
                                readmode,
                                prefixlength,
                                numofparts,
                                maximumspace,
                                NULL,
                                NULL,
                                sfxstrategy,
                                sfxprogress,
                                withprogressbar,
                                logger,
                                err);
}

static void gt_suffixer_sort_with_dcov(void *voidsfi,
                                       GtSuffixsortspace *sssp,
                                       GtUword blisbl,
                                       GtUword width,
                                       GtUword depth)
{
  Sfxiterator *sfi = (Sfxiterator *) voidsfi;

  gt_assert(sfi != NULL);
  gt_differencecover_sortunsortedbucket(sssp,
                                        gt_Outlcpinfo_lcpvalues_ref(
                                               sfi->outlcpinfo),
                                        sfi->dcov, blisbl, width,depth);
}

static void gt_sfxiterator_preparethispart(Sfxiterator *sfi)
{
  GtUword sumofwidthforpart;
  GtBucketspec2 *bucketspec2 = NULL;

  if (sfi->part == 0 && sfi->withprogressbar)
  {
    gt_assert(sfi->bcktab != NULL);
    gt_progressbar_start(&sfi->bucketiterstep,
                         (GtUint64)
                         gt_bcktab_numofallcodes(sfi->bcktab));
  }
  sfi->currentmincode = gt_suftabparts_minindex(sfi->part,sfi->suftabparts);
  sfi->currentmaxcode = gt_suftabparts_maxindex(sfi->part,sfi->suftabparts);
  sfi->widthofpart = gt_suftabparts_widthofpart(sfi->part,sfi->suftabparts);
  if (sfi->sfxprogress != NULL)
  {
    gt_timer_show_progress(sfi->sfxprogress, "inserting suffixes into buckets",
                           stdout);
  }
  if (gt_suftabparts_numofparts(sfi->suftabparts) > 1U)
  {
    gt_logger_log(sfi->logger,"compute part %u: "
                              ""GT_WU" suffixes,"GT_WU" buckets from "
                              ""GT_WU".."GT_WU"",
                  sfi->part,
                  gt_suftabparts_widthofpart(sfi->part,sfi->suftabparts),
                  sfi->currentmaxcode - sfi->currentmincode + 1,
                  sfi->currentmincode,
                  sfi->currentmaxcode);
    gt_bcktab_assignboundsforpart(sfi->bcktab,
                                  sfi->currentmincode,
                                  sfi->currentmaxcode);
    if (sfi->mappedmarkprefixbuckets != NULL)
    {
      sfi->markprefixbuckets
        = (GtBitsequence *)
          gt_Sfxmappedrange_map(sfi->mappedmarkprefixbuckets,
                                sfi->currentmincode,
                                sfi->currentmaxcode);
    }
  }
  SHOWACTUALSPACE;
  gt_suffixsortspace_partoffset_set(sfi->suffixsortspace,
                                    gt_suftabparts_offset(sfi->part,
                                                          sfi->suftabparts));
  if (sfi->sfxstrategy.spmopt_minlength == 0)
  {
    if (sfi->sfxstrategy.storespecialcodes)
    {
      sfx_derivespecialcodesfromtable(sfi,
                          gt_suftabparts_numofparts(sfi->suftabparts) == 1U
                            ? true
                            : false);
    } else
    {
      sfx_derivespecialcodesonthefly(sfi);
    }
  }
  SHOWACTUALSPACE;
  sfi->exportptr = gt_suffixsortspace_exportptr(sfi->suffixsortspace, 0);
  if (sfi->prefixlength > 1U
      && gt_encseq_has_twobitencoding(sfi->encseq)
      && !sfi->sfxstrategy.kmerswithencseqreader)
  {
    insertsuffix_getencseqkmers_twobitencoding(
                                     sfi->encseq,
                                     sfi->readmode,
                                     sfi->sfxstrategy.spmopt_minlength == 0
                                       ? sfi->prefixlength
                                       : sfi->spmopt_kmerscansize,
                                     sfi->sfxstrategy.spmopt_minlength == 0
                                       ? sfi->prefixlength
                                       : sfi->sfxstrategy.spmopt_minlength,
                                     sfi,
                                     NULL);
  } else
  {
    if (sfi->sfxstrategy.iteratorbasedkmerscanning)
    {
      getencseqkmersinsertkmerwithoutspecial(sfi->encseq,
                                             sfi->readmode,
                                             sfi->prefixlength,
                                             sfi);
    } else
    {
      getencseqkmers(sfi->encseq,sfi->readmode,sfi->prefixlength,
                     gt_insertkmerwithoutspecial,sfi);
    }
  }
  SHOWACTUALSPACE;
  gt_suffixsortspace_export_done(sfi->suffixsortspace);
  if (sfi->sfxprogress != NULL)
  {
    gt_timer_show_progress(sfi->sfxprogress, "sorting the buckets", stdout);
  }
  /* exit(0); just for testing */
  sumofwidthforpart = gt_suftabparts_sumofwidth(sfi->part,sfi->suftabparts);

  if (gt_suftabparts_numofparts(sfi->suftabparts) == 1U &&
      sfi->outlcpinfo == NULL &&
      sfi->prefixlength >= 2U &&
      sfi->sfxstrategy.spmopt_minlength == 0 &&
      GT_SFX_THREADS_JOBS == 1U)
  {
    bucketspec2 = gt_copysort_new(sfi->bcktab,sfi->encseq,sfi->readmode,
                                  sumofwidthforpart,sfi->numofchars);
  }
  SHOWACTUALSPACE;
  if (sfi->part == 0)
  {
    gt_logger_log(sfi->logger,"used workspace for sorting: %.2f MB",
                  GT_MEGABYTES(gt_size_of_sort_workspace (&sfi->sfxstrategy)));
  }
  if (!sfi->sfxstrategy.onlybucketinsertion)
  {
    unsigned int sortmaxdepth;
    GtProcessunsortedsuffixrange processunsortedsuffixrange;
    void *processunsortedsuffixrangeinfo;

    if (sfi->dcov == NULL)
    {
      if (sfi->sfxstrategy.userdefinedsortmaxdepth == 0)
      {
        sortmaxdepth = 0;
      } else
      {
        sortmaxdepth = sfi->sfxstrategy.userdefinedsortmaxdepth;
      }
      processunsortedsuffixrange = NULL;
      processunsortedsuffixrangeinfo = NULL;
    } else
    {
      gt_assert(sfi->sfxstrategy.userdefinedsortmaxdepth == 0);
      sortmaxdepth = sfi->sfxstrategy.differencecover;
      processunsortedsuffixrange = gt_suffixer_sort_with_dcov;
      processunsortedsuffixrangeinfo = sfi;
    }
    gt_assert(sortmaxdepth != 0 || processunsortedsuffixrange == NULL);
    gt_bcktab_determinemaxsize(sfi->bcktab, sfi->currentmincode,
                               sfi->currentmaxcode,sumofwidthforpart);
#ifdef GT_THREADS_ENABLED
    if (GT_SFX_THREADS_JOBS > 1U
#ifdef GT_THREADS_PARTITION
        &&
        sfi->partitions_for_threads != NULL &&
        gt_suftabparts_numofparts(sfi->partitions_for_threads[sfi->part]) > 1U
#endif
        )
    {
#ifdef GT_THREADS_PARTITION
      gt_threaded_partition_sortallbuckets
                                 (sfi->suffixsortspace,
                                 sfi->partitions_for_threads[sfi->part],
                                 sfi->encseq,
                                 sfi->readmode,
                                 sfi->bcktab,
                                 sfi->numofchars,
                                 sfi->prefixlength,
                                 sortmaxdepth,
                                 &sfi->sfxstrategy,
                                 processunsortedsuffixrange,
                                 processunsortedsuffixrangeinfo,
                                 sfi->logger);
#else
      gt_threaded_stream_sortallbuckets
                                 (sfi->suffixsortspace,
                                 sfi->encseq,
                                 sfi->readmode,
                                 sfi->bcktab,
                                 sfi->currentmincode,
                                 sfi->currentmaxcode,
                                 sumofwidthforpart,
                                 sfi->numofchars,
                                 sfi->prefixlength,
                                 sortmaxdepth,
                                 &sfi->sfxstrategy,
                                 processunsortedsuffixrange,
                                 processunsortedsuffixrangeinfo,
                                 sfi->logger);
#endif
    } else
    {
#endif
    gt_sortallbuckets(sfi->suffixsortspace,
                      sumofwidthforpart,
                      bucketspec2,
                      sfi->encseq,
                      sfi->readmode,
                      sfi->currentmincode,
                      sfi->currentmaxcode,
                      sfi->bcktab,
                      sfi->numofchars,
                      sfi->prefixlength,
                      sfi->outlcpinfo,
                      sortmaxdepth,
                      &sfi->sfxstrategy,
                      processunsortedsuffixrange,
                      processunsortedsuffixrangeinfo,
                      &sfi->bucketiterstep,
                      sfi->logger);
#ifdef GT_THREADS_ENABLED
    }
#endif
  }
  if (bucketspec2 != NULL)
  {
    gt_copysort_derivesorting(bucketspec2,sfi->suffixsortspace,sfi->logger);
    gt_copysort_delete(bucketspec2);
    bucketspec2 = NULL;
  }
  SHOWACTUALSPACE;
  sfi->part++;
}

const GtSuffixsortspace *gt_Sfxiterator_next(GtUword *numberofsuffixes,
                                             bool *specialsuffixes,
                                             Sfxiterator *sfi)
{
  if (sfi->part < gt_suftabparts_numofparts(sfi->suftabparts))
  {
    gt_sfxiterator_preparethispart(sfi);
    *numberofsuffixes = sfi->widthofpart;
    if (specialsuffixes != NULL)
    {
      *specialsuffixes = false;
    }
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
  if (sfi->sfxstrategy.spmopt_minlength == 0)
  {
    gt_suffixsortspace_partoffset_set(sfi->suffixsortspace,
                                      sfi->totallength-sfi->specialcharacters);
    sfi->exhausted = gt_SSSPbuf_fillspecialnextpage(sfi->suffixsortspace,
                                                    sfi->readmode,
                                                    sfi->sri,
                                                    sfi->totallength,
                                                    sfi->sssp_buf);
    *numberofsuffixes = gt_SSSPbuf_filled(sfi->sssp_buf);
  } else
  {
    sfi->exhausted = true;
    *numberofsuffixes = 0;
  }
  if (specialsuffixes != NULL)
  {
    *specialsuffixes = true;
  }
  return sfi->suffixsortspace;
}

int gt_Sfxiterator_bcktab2file(FILE *fp, Sfxiterator *sfi, GtError *err)
{
  gt_error_check(err);
  gt_assert(sfi != NULL && sfi->bcktab != NULL);
  if (gt_suftabparts_numofparts(sfi->suftabparts) <= 1U)
  {
    int ret = gt_bcktab_flush_to_file(fp,sfi->bcktab,err);
    gt_fa_fclose(fp);
    return ret;
  }
  return 0;
}

GtUword gt_Sfxiterator_longest(const Sfxiterator *sfi)
{
  gt_assert(sfi != NULL);
  if (sfi->sfxstrategy.spmopt_minlength == 0)
  {
    return gt_suffixsortspace_longest(sfi->suffixsortspace);
  } else
  {
    return 0;
  }
}
