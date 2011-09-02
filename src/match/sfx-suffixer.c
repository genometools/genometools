/*
  Copyright (c) 2007-2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2011 Center for Bioinformatics, University of Hamburg

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
#include "extended/uint64hashtable.h"
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
#include "sfx-maprange.h"
#include "stamp.h"

typedef struct
{
  unsigned long allocatedSuffixptr, nextfreeSuffixptr;
  GtSuffixsortspace *sssp;
} GtSuffixposbuffer;

struct Sfxiterator
{
  /* globally constant */
  const GtEncseq *encseq;
  GtReadmode readmode;
  unsigned long specialcharacters,
                totallength;
  unsigned int numofchars,
               prefixlength;
  Sfxstrategy sfxstrategy;
  bool withprogressbar;

  /* invariant for each part */
  Suftabparts *suftabparts;
  Outlcpinfo *outlcpinfoforsample;
  GtBcktab *bcktab;
  GtLeftborder *leftborder; /* points to bcktab->leftborder */
  Differencecover *dcov;

  /* changed in each part */
  GtSuffixsortspace *suffixsortspace;
  GtCodetype currentmincode,
             currentmaxcode;
  unsigned long widthofpart;
  unsigned int part;
  Outlcpinfo *outlcpinfo;
  GtSuffixposbuffer fusp;
  GtRange overhang;
  bool exhausted;
  unsigned long long bucketiterstep; /* for progressbar */
  GtLogger *logger;
  GtTimer *sfxprogress;
  GtSpecialrangeiterator *sri; /* refers to space used in each part */

  /* use for generating k-mer codes */
  FILE *outfpbcktab;
  bool storespecials;
  unsigned long nextfreeCodeatposition;
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
  unsigned long spmopt_numofallprefixcodes,
                spmopt_numofallsuffixcodes;
  GtSfxmappedrange *mappedmarkprefixbuckets;
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

static void compareCodeatpositionlists(const Codeatposition *codelist1,
                                       unsigned long len1,
                                       const Codeatposition *codelist2,
                                       unsigned long len2)
{
  unsigned long idx;

  if (len1 != len2)
  {
    fprintf(stderr,"%s: len1 = %lu != %lu = len2\n",__func__,len1,len2);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  for (idx=0; idx<len1; idx++)
  {
    if (codelist1[idx].position != codelist2[idx].position)
    {
      fprintf(stderr,"%s: listlength = %lu,idx %lu, codelist1.position "
                     "= %lu != %lu = codelist2.position\n",
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
      fprintf(stderr,"%s: idx %lu, codelist1.maxprefixindex = %u != %u = "
                     "codelist2.maxprefixindex\n",__func__,idx,
                      codelist1[idx].maxprefixindex,
                      codelist2[idx].maxprefixindex);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    if (codelist1[idx].code != codelist2[idx].code)
    {
      fprintf(stderr,"%s: idx %lu, codelist1.code = %u != %u = "
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
  gt_assert(realspecialranges+1 >= nextfreeCodeatposition2);
  compareCodeatpositionlists(spaceCodeatposition1,
                             nextfreeCodeatposition1,
                             spaceCodeatposition2,
                             nextfreeCodeatposition2);
  gt_free(spaceCodeatposition2);
}

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
        fprintf(stderr,"%s: ### position %lu, code2 = %lu != 0\n",__func__,
                       position,code2);
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

#define GT_SCANCODE_TO_BCKCODE(CODE)\
        (GtCodetype) ((CODE) >> sfi->spmopt_kmerscancodeshift2bckcode)

#define GT_SCANCODE_TO_PREFIXCODE(CODE)\
        (GtCodetype) ((CODE) >> sfi->spmopt_kmerscancodeshift2prefixcode)

#ifdef _LP64
#define GT_SCANCODE_TO_SUFFIXCODE(CODE)\
        (GtCodetype) ((CODE) & sfi->spmopt_kmerscancodesuffixmask);

static bool gt_checksuffixprefixbuckets(const Sfxiterator *sfi,
                                        GtCodetype scancode)
{
  GtCodetype prefixcode = GT_SCANCODE_TO_PREFIXCODE(scancode);
  GtCodetype suffixcode = GT_SCANCODE_TO_SUFFIXCODE(scancode);

  gt_assert(prefixcode < sfi->spmopt_numofallprefixcodes);
  gt_assert(suffixcode < sfi->spmopt_numofallsuffixcodes);
  return (GT_ISIBITSET(sfi->markprefixbuckets,prefixcode) &&
          GT_ISIBITSET(sfi->marksuffixbuckets,suffixcode)) ? true : false;
}
#else
static bool gt_checksuffixprefixbuckets(const Sfxiterator *sfi,
                                        GtCodetype scancode)
{
  GtCodetype prefixcode = GT_SCANCODE_TO_PREFIXCODE(scancode);

  gt_assert(prefixcode < sfi->spmopt_numofallprefixcodes);
  return GT_ISIBITSET(sfi->markprefixbuckets,prefixcode) ? true : false;
}
#endif

static void gt_insertkmerwithoutspecial1(void *processinfo,
                                         bool firstinrange,
                                         unsigned long position,
                                         GtCodetype scancode)
{
  Sfxiterator *sfi = (Sfxiterator *) processinfo;

  if (sfi->markprefixbuckets == NULL)
  {
    if (scancode >= sfi->currentmincode && scancode <= sfi->currentmaxcode)
    {
      unsigned long stidx;

      stidx = gt_bcktab_leftborder_insertionindex(sfi->leftborder,scancode);
      /* from right to left */
      GT_SUFFIXSORTSPACE_EXPORT_SET(sfi->suffixsortspace,sfi->exportptr,stidx,
                                    position);
    }
  } else
  {
    GtCodetype bcktabcode = GT_SCANCODE_TO_BCKCODE(scancode);

    if (bcktabcode >= sfi->currentmincode &&
        bcktabcode <= sfi->currentmaxcode &&
        (firstinrange || gt_checksuffixprefixbuckets(sfi,scancode)))
    {
      unsigned long stidx;

      stidx = gt_bcktab_leftborder_insertionindex(sfi->leftborder,bcktabcode);
      /* from right to left */
      GT_SUFFIXSORTSPACE_EXPORT_SET(sfi->suffixsortspace,sfi->exportptr,stidx,
                                    position);
    }
  }
}

static void gt_insertkmerwithoutspecial(void *processinfo,
                                        unsigned long position,
                                        const GtKmercode *kmercode)
{
  if (!kmercode->definedspecialposition)
  {
    gt_insertkmerwithoutspecial1(processinfo,false, position, kmercode->code);
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
  if (sfi->suftabparts != NULL && stpgetnumofparts(sfi->suftabparts) > 1U &&
      sfi->outfpbcktab != NULL)
  {
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
  gt_bcktab_delete(sfi->bcktab,sfi->logger);
  gt_freesuftabparts(sfi->suftabparts);
  gt_Outlcpinfo_delete(sfi->outlcpinfoforsample);
  if (sfi->mappedmarkprefixbuckets == NULL)
  {
    gt_free(sfi->markprefixbuckets);
  }
  gt_Sfxmappedrange_delete(sfi->mappedmarkprefixbuckets,sfi->logger);
  sfi->mappedmarkprefixbuckets = NULL;
  gt_free(sfi->marksuffixbuckets);
  gt_differencecover_delete(sfi->dcov);
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
    unsigned long position = 0;

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
    unsigned long position = 0;

    while ((kmercodeptr = gt_kmercodeiterator_encseq_next(kmercodeiterator))
                          != NULL)
    {
      gt_insertkmerwithoutspecial(sfi,position++,kmercodeptr);
    }
    gt_kmercodeiterator_delete(kmercodeiterator);
  }
}

static int computepartsfittingmaximumspace(size_t estimatedspace,
                                           unsigned long maximumspace,
                                           const GtBcktab *bcktab,
                                           const GtSfxmappedrange
                                             *mappedmarkprefixbuckets,
                                           unsigned long totallength,
                                           unsigned long specialcharacters,
                                           unsigned long numofsuffixestosort,
                                           bool suftabuint,
                                           GtError *err)
{
  unsigned int parts;
  Suftabparts *suftabparts;
  unsigned long size_lb_cs = gt_bcktab_size_lb_cs(bcktab);
  size_t sizeofprefixmarks;

  gt_error_check(err);
  if (mappedmarkprefixbuckets != NULL)
  {
    sizeofprefixmarks = gt_Sfxmappedrange_size_entire(mappedmarkprefixbuckets);
  } else
  {
    sizeofprefixmarks = 0;
  }
  /*
  printf("maxspace=%.2f\n",GT_MEGABYTES(maximumspace));
  printf("estimatedspace=%.2f\n",GT_MEGABYTES(estimatedspace));
  */
  if (estimatedspace >= (size_t) maximumspace)
  {
    gt_error_set(err,"already used %.2f MB of memory, cannot compute "
                     "index in at most %.2f MB",
                     GT_MEGABYTES(estimatedspace), GT_MEGABYTES(maximumspace));
    return -1;
  }
  for (parts = 1U; parts <= 500U; parts++)
  {
    size_t suftabsize;

    suftabparts = gt_newsuftabparts(parts,
                                    bcktab,
                                    mappedmarkprefixbuckets,
                                    numofsuffixestosort,
                                    specialcharacters + 1,
                                    NULL);
    gt_assert(suftabparts != NULL);
    suftabsize = gt_suffixsortspace_requiredspace(
                                         stpgetlargestsuftabwidth(suftabparts),
                                         totallength,
                                         suftabuint);
    /*
    printf("parts=%u,suftabsize=%.2f,largestsizemappedpartwise=%.2f,"
           "size_lb_cs=%.2f,sizeprefixmarks=%.2f\n",
           parts,
           GT_MEGABYTES(suftabsize),
           GT_MEGABYTES(stpgetlargestsizemappedpartwise(suftabparts)),
           GT_MEGABYTES(size_lb_cs),
           GT_MEGABYTES(sizeofprefixmarks));
    */
    if (parts == 1U)
    {
      if ((unsigned long) (suftabsize + estimatedspace) <= maximumspace)
      {
        gt_freesuftabparts(suftabparts);
        return (int) parts;
      }
    } else
    {
      if ((unsigned long) (suftabsize +
                           stpgetlargestsizemappedpartwise(suftabparts) +
                           estimatedspace - size_lb_cs - sizeofprefixmarks)
                           <= maximumspace)
      {
        gt_freesuftabparts(suftabparts);
        return (int) parts;
      }
    }
    gt_freesuftabparts(suftabparts);
  }
  gt_error_set(err,"cannot compute enhanced suffix array in at most %lu bytes",
                   maximumspace);
  return -1;
}

#undef DEBUGSIZEESTIMATION
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

GtCodetype gt_kmercode_at_position(const GtTwobitencoding *twobitencoding,
                                   unsigned long pos,
                                   unsigned int kmersize)
{
  const unsigned int unitoffset = (unsigned int) GT_MODBYUNITSIN2BITENC(pos);
  const unsigned long unitindex = GT_DIVBYUNITSIN2BITENC(pos);
  const GtCodetype maskright = GT_MASKRIGHT(kmersize);

  if (unitoffset <= (unsigned int) GT_UNITSIN2BITENC - kmersize)
  {
    return (GtCodetype) (twobitencoding[unitindex]
            >> GT_MULT2(GT_UNITSIN2BITENC - kmersize - unitoffset))
           & maskright;
  } else
  {
    unsigned int shiftleft = GT_MULT2(unitoffset+kmersize-GT_UNITSIN2BITENC);
    return (GtCodetype)
           ((twobitencoding[unitindex] << shiftleft) |
            (twobitencoding[unitindex+1] >> (GT_MULT2(GT_UNITSIN2BITENC) -
                                             shiftleft)))
           & maskright;
  }
}

GtCodetype gt_kmercode_at_firstpos(const GtTwobitencoding *twobitencoding,
                                   unsigned int kmersize)
{
  const GtCodetype maskright = GT_MASKRIGHT(kmersize);
  return (GtCodetype) (twobitencoding[0] >>
                       GT_MULT2(GT_UNITSIN2BITENC - kmersize)) & maskright;
}

#define GT_SWAPBITPAIRS(KMER,L1,L2,D) (((KMER) & (3UL << L1)) >> D) |\
                                      (((KMER) & (3UL << L2)) << D)

static GtCodetype gt_kmercode_reverse(GtCodetype kmer,unsigned int kmersize)
{
  switch (kmersize)
  {
    case 2:
      return GT_SWAPBITPAIRS(kmer,2,0,2);
    case 3:
      return GT_SWAPBITPAIRS(kmer,4,0,4) |
             (kmer & (3U << 2));
    case 4:
      return GT_SWAPBITPAIRS(kmer,6,0,6) |
             GT_SWAPBITPAIRS(kmer,4,2,2);
    case 5:
      return GT_SWAPBITPAIRS(kmer,8,0,8) |
             GT_SWAPBITPAIRS(kmer,6,2,4) |
             (kmer & (3U << 4));
    case 6:
      return GT_SWAPBITPAIRS(kmer,10,0,10) |
             GT_SWAPBITPAIRS(kmer,8,2,6) |
             GT_SWAPBITPAIRS(kmer,6,4,2);
    case 7:
      return GT_SWAPBITPAIRS(kmer,12,0,12) |
             GT_SWAPBITPAIRS(kmer,10,2,8) |
             GT_SWAPBITPAIRS(kmer,8,4,4) |
             (kmer & (3U << 6));
    case 8:
      return GT_SWAPBITPAIRS(kmer,14,0,14) |
             GT_SWAPBITPAIRS(kmer,12,2,10) |
             GT_SWAPBITPAIRS(kmer,10,4,6) |
             GT_SWAPBITPAIRS(kmer,8,6,2);
    case 9:
      return GT_SWAPBITPAIRS(kmer,16,0,16) |
             GT_SWAPBITPAIRS(kmer,14,2,12) |
             GT_SWAPBITPAIRS(kmer,12,4,8) |
             GT_SWAPBITPAIRS(kmer,10,6,4) |
             (kmer & (3U << 8));
    case 10:
      return GT_SWAPBITPAIRS(kmer,18,0,18) |
             GT_SWAPBITPAIRS(kmer,16,2,14) |
             GT_SWAPBITPAIRS(kmer,14,4,10) |
             GT_SWAPBITPAIRS(kmer,12,6,6) |
             GT_SWAPBITPAIRS(kmer,10,8,2);
    case 11:
      return GT_SWAPBITPAIRS(kmer,20,0,20) |
             GT_SWAPBITPAIRS(kmer,18,2,16) |
             GT_SWAPBITPAIRS(kmer,16,4,12) |
             GT_SWAPBITPAIRS(kmer,14,6,8) |
             GT_SWAPBITPAIRS(kmer,12,8,4) |
             (kmer & (3U << 10));
    case 12:
      return GT_SWAPBITPAIRS(kmer,22,0,22) |
             GT_SWAPBITPAIRS(kmer,20,2,18) |
             GT_SWAPBITPAIRS(kmer,18,4,14) |
             GT_SWAPBITPAIRS(kmer,16,6,10) |
             GT_SWAPBITPAIRS(kmer,14,8,6) |
             GT_SWAPBITPAIRS(kmer,12,10,2);
    case 13:
      return GT_SWAPBITPAIRS(kmer,24,0,24) |
             GT_SWAPBITPAIRS(kmer,22,2,20) |
             GT_SWAPBITPAIRS(kmer,20,4,16) |
             GT_SWAPBITPAIRS(kmer,18,6,12) |
             GT_SWAPBITPAIRS(kmer,16,8,8) |
             GT_SWAPBITPAIRS(kmer,14,10,4) |
             (kmer & (3U << 12));
    case 14:
      return GT_SWAPBITPAIRS(kmer,26,0,26) |
             GT_SWAPBITPAIRS(kmer,24,2,22) |
             GT_SWAPBITPAIRS(kmer,22,4,18) |
             GT_SWAPBITPAIRS(kmer,20,6,14) |
             GT_SWAPBITPAIRS(kmer,18,8,10) |
             GT_SWAPBITPAIRS(kmer,16,10,6) |
             GT_SWAPBITPAIRS(kmer,14,12,2);
    case 15:
      return GT_SWAPBITPAIRS(kmer,28,0,28) |
             GT_SWAPBITPAIRS(kmer,26,2,24) |
             GT_SWAPBITPAIRS(kmer,24,4,20) |
             GT_SWAPBITPAIRS(kmer,22,6,16) |
             GT_SWAPBITPAIRS(kmer,20,8,12) |
             GT_SWAPBITPAIRS(kmer,18,10,8) |
             GT_SWAPBITPAIRS(kmer,16,12,4) |
             (kmer & (3U << 14));
#ifdef _LP64
    case 16:
      return GT_SWAPBITPAIRS(kmer,30,0,30) |
             GT_SWAPBITPAIRS(kmer,28,2,26) |
             GT_SWAPBITPAIRS(kmer,26,4,22) |
             GT_SWAPBITPAIRS(kmer,24,6,18) |
             GT_SWAPBITPAIRS(kmer,22,8,14) |
             GT_SWAPBITPAIRS(kmer,20,10,10) |
             GT_SWAPBITPAIRS(kmer,18,12,6) |
             GT_SWAPBITPAIRS(kmer,16,14,2);
    case 17:
      return GT_SWAPBITPAIRS(kmer,32,0,32) |
             GT_SWAPBITPAIRS(kmer,30,2,28) |
             GT_SWAPBITPAIRS(kmer,28,4,24) |
             GT_SWAPBITPAIRS(kmer,26,6,20) |
             GT_SWAPBITPAIRS(kmer,24,8,16) |
             GT_SWAPBITPAIRS(kmer,22,10,12) |
             GT_SWAPBITPAIRS(kmer,20,12,8) |
             GT_SWAPBITPAIRS(kmer,18,14,4) |
             (kmer & (3U << 16));
    case 18:
      return GT_SWAPBITPAIRS(kmer,34,0,34) |
             GT_SWAPBITPAIRS(kmer,32,2,30) |
             GT_SWAPBITPAIRS(kmer,30,4,26) |
             GT_SWAPBITPAIRS(kmer,28,6,22) |
             GT_SWAPBITPAIRS(kmer,26,8,18) |
             GT_SWAPBITPAIRS(kmer,24,10,14) |
             GT_SWAPBITPAIRS(kmer,22,12,10) |
             GT_SWAPBITPAIRS(kmer,20,14,6) |
             GT_SWAPBITPAIRS(kmer,18,16,2);
    case 19:
      return GT_SWAPBITPAIRS(kmer,36,0,36) |
             GT_SWAPBITPAIRS(kmer,34,2,32) |
             GT_SWAPBITPAIRS(kmer,32,4,28) |
             GT_SWAPBITPAIRS(kmer,30,6,24) |
             GT_SWAPBITPAIRS(kmer,28,8,20) |
             GT_SWAPBITPAIRS(kmer,26,10,16) |
             GT_SWAPBITPAIRS(kmer,24,12,12) |
             GT_SWAPBITPAIRS(kmer,22,14,8) |
             GT_SWAPBITPAIRS(kmer,20,16,4) |
             (kmer & (3U << 18));
    case 20:
      return GT_SWAPBITPAIRS(kmer,38,0,38) |
             GT_SWAPBITPAIRS(kmer,36,2,34) |
             GT_SWAPBITPAIRS(kmer,34,4,30) |
             GT_SWAPBITPAIRS(kmer,32,6,26) |
             GT_SWAPBITPAIRS(kmer,30,8,22) |
             GT_SWAPBITPAIRS(kmer,28,10,18) |
             GT_SWAPBITPAIRS(kmer,26,12,14) |
             GT_SWAPBITPAIRS(kmer,24,14,10) |
             GT_SWAPBITPAIRS(kmer,22,16,6) |
             GT_SWAPBITPAIRS(kmer,20,18,2);
    case 21:
      return GT_SWAPBITPAIRS(kmer,40,0,40) |
             GT_SWAPBITPAIRS(kmer,38,2,36) |
             GT_SWAPBITPAIRS(kmer,36,4,32) |
             GT_SWAPBITPAIRS(kmer,34,6,28) |
             GT_SWAPBITPAIRS(kmer,32,8,24) |
             GT_SWAPBITPAIRS(kmer,30,10,20) |
             GT_SWAPBITPAIRS(kmer,28,12,16) |
             GT_SWAPBITPAIRS(kmer,26,14,12) |
             GT_SWAPBITPAIRS(kmer,24,16,8) |
             GT_SWAPBITPAIRS(kmer,22,18,4) |
             (kmer & (3U << 20));
    case 22:
      return GT_SWAPBITPAIRS(kmer,42,0,42) |
             GT_SWAPBITPAIRS(kmer,40,2,38) |
             GT_SWAPBITPAIRS(kmer,38,4,34) |
             GT_SWAPBITPAIRS(kmer,36,6,30) |
             GT_SWAPBITPAIRS(kmer,34,8,26) |
             GT_SWAPBITPAIRS(kmer,32,10,22) |
             GT_SWAPBITPAIRS(kmer,30,12,18) |
             GT_SWAPBITPAIRS(kmer,28,14,14) |
             GT_SWAPBITPAIRS(kmer,26,16,10) |
             GT_SWAPBITPAIRS(kmer,24,18,6) |
             GT_SWAPBITPAIRS(kmer,22,20,2);
    case 23:
      return GT_SWAPBITPAIRS(kmer,44,0,44) |
             GT_SWAPBITPAIRS(kmer,42,2,40) |
             GT_SWAPBITPAIRS(kmer,40,4,36) |
             GT_SWAPBITPAIRS(kmer,38,6,32) |
             GT_SWAPBITPAIRS(kmer,36,8,28) |
             GT_SWAPBITPAIRS(kmer,34,10,24) |
             GT_SWAPBITPAIRS(kmer,32,12,20) |
             GT_SWAPBITPAIRS(kmer,30,14,16) |
             GT_SWAPBITPAIRS(kmer,28,16,12) |
             GT_SWAPBITPAIRS(kmer,26,18,8) |
             GT_SWAPBITPAIRS(kmer,24,20,4) |
             (kmer & (3U << 22));
    case 24:
      return GT_SWAPBITPAIRS(kmer,46,0,46) |
             GT_SWAPBITPAIRS(kmer,44,2,42) |
             GT_SWAPBITPAIRS(kmer,42,4,38) |
             GT_SWAPBITPAIRS(kmer,40,6,34) |
             GT_SWAPBITPAIRS(kmer,38,8,30) |
             GT_SWAPBITPAIRS(kmer,36,10,26) |
             GT_SWAPBITPAIRS(kmer,34,12,22) |
             GT_SWAPBITPAIRS(kmer,32,14,18) |
             GT_SWAPBITPAIRS(kmer,30,16,14) |
             GT_SWAPBITPAIRS(kmer,28,18,10) |
             GT_SWAPBITPAIRS(kmer,26,20,6) |
             GT_SWAPBITPAIRS(kmer,24,22,2);
    case 25:
      return GT_SWAPBITPAIRS(kmer,48,0,48) |
             GT_SWAPBITPAIRS(kmer,46,2,44) |
             GT_SWAPBITPAIRS(kmer,44,4,40) |
             GT_SWAPBITPAIRS(kmer,42,6,36) |
             GT_SWAPBITPAIRS(kmer,40,8,32) |
             GT_SWAPBITPAIRS(kmer,38,10,28) |
             GT_SWAPBITPAIRS(kmer,36,12,24) |
             GT_SWAPBITPAIRS(kmer,34,14,20) |
             GT_SWAPBITPAIRS(kmer,32,16,16) |
             GT_SWAPBITPAIRS(kmer,30,18,12) |
             GT_SWAPBITPAIRS(kmer,28,20,8) |
             GT_SWAPBITPAIRS(kmer,26,22,4) |
             (kmer & (3U << 24));
    case 26:
      return GT_SWAPBITPAIRS(kmer,50,0,50) |
             GT_SWAPBITPAIRS(kmer,48,2,46) |
             GT_SWAPBITPAIRS(kmer,46,4,42) |
             GT_SWAPBITPAIRS(kmer,44,6,38) |
             GT_SWAPBITPAIRS(kmer,42,8,34) |
             GT_SWAPBITPAIRS(kmer,40,10,30) |
             GT_SWAPBITPAIRS(kmer,38,12,26) |
             GT_SWAPBITPAIRS(kmer,36,14,22) |
             GT_SWAPBITPAIRS(kmer,34,16,18) |
             GT_SWAPBITPAIRS(kmer,32,18,14) |
             GT_SWAPBITPAIRS(kmer,30,20,10) |
             GT_SWAPBITPAIRS(kmer,28,22,6) |
             GT_SWAPBITPAIRS(kmer,26,24,2);
    case 27:
      return GT_SWAPBITPAIRS(kmer,52,0,52) |
             GT_SWAPBITPAIRS(kmer,50,2,48) |
             GT_SWAPBITPAIRS(kmer,48,4,44) |
             GT_SWAPBITPAIRS(kmer,46,6,40) |
             GT_SWAPBITPAIRS(kmer,44,8,36) |
             GT_SWAPBITPAIRS(kmer,42,10,32) |
             GT_SWAPBITPAIRS(kmer,40,12,28) |
             GT_SWAPBITPAIRS(kmer,38,14,24) |
             GT_SWAPBITPAIRS(kmer,36,16,20) |
             GT_SWAPBITPAIRS(kmer,34,18,16) |
             GT_SWAPBITPAIRS(kmer,32,20,12) |
             GT_SWAPBITPAIRS(kmer,30,22,8) |
             GT_SWAPBITPAIRS(kmer,28,24,4) |
             (kmer & (3U << 26));
    case 28:
      return GT_SWAPBITPAIRS(kmer,54,0,54) |
             GT_SWAPBITPAIRS(kmer,52,2,50) |
             GT_SWAPBITPAIRS(kmer,50,4,46) |
             GT_SWAPBITPAIRS(kmer,48,6,42) |
             GT_SWAPBITPAIRS(kmer,46,8,38) |
             GT_SWAPBITPAIRS(kmer,44,10,34) |
             GT_SWAPBITPAIRS(kmer,42,12,30) |
             GT_SWAPBITPAIRS(kmer,40,14,26) |
             GT_SWAPBITPAIRS(kmer,38,16,22) |
             GT_SWAPBITPAIRS(kmer,36,18,18) |
             GT_SWAPBITPAIRS(kmer,34,20,14) |
             GT_SWAPBITPAIRS(kmer,32,22,10) |
             GT_SWAPBITPAIRS(kmer,30,24,6) |
             GT_SWAPBITPAIRS(kmer,28,26,2);
    case 29:
      return GT_SWAPBITPAIRS(kmer,56,0,56) |
             GT_SWAPBITPAIRS(kmer,54,2,52) |
             GT_SWAPBITPAIRS(kmer,52,4,48) |
             GT_SWAPBITPAIRS(kmer,50,6,44) |
             GT_SWAPBITPAIRS(kmer,48,8,40) |
             GT_SWAPBITPAIRS(kmer,46,10,36) |
             GT_SWAPBITPAIRS(kmer,44,12,32) |
             GT_SWAPBITPAIRS(kmer,42,14,28) |
             GT_SWAPBITPAIRS(kmer,40,16,24) |
             GT_SWAPBITPAIRS(kmer,38,18,20) |
             GT_SWAPBITPAIRS(kmer,36,20,16) |
             GT_SWAPBITPAIRS(kmer,34,22,12) |
             GT_SWAPBITPAIRS(kmer,32,24,8) |
             GT_SWAPBITPAIRS(kmer,30,26,4) |
             (kmer & (3U << 28));
    case 30:
      return GT_SWAPBITPAIRS(kmer,58,0,58) |
             GT_SWAPBITPAIRS(kmer,56,2,54) |
             GT_SWAPBITPAIRS(kmer,54,4,50) |
             GT_SWAPBITPAIRS(kmer,52,6,46) |
             GT_SWAPBITPAIRS(kmer,50,8,42) |
             GT_SWAPBITPAIRS(kmer,48,10,38) |
             GT_SWAPBITPAIRS(kmer,46,12,34) |
             GT_SWAPBITPAIRS(kmer,44,14,30) |
             GT_SWAPBITPAIRS(kmer,42,16,26) |
             GT_SWAPBITPAIRS(kmer,40,18,22) |
             GT_SWAPBITPAIRS(kmer,38,20,18) |
             GT_SWAPBITPAIRS(kmer,36,22,14) |
             GT_SWAPBITPAIRS(kmer,34,24,10) |
             GT_SWAPBITPAIRS(kmer,32,26,6) |
             GT_SWAPBITPAIRS(kmer,30,28,2);
    case 31:
      return GT_SWAPBITPAIRS(kmer,60,0,60) |
             GT_SWAPBITPAIRS(kmer,58,2,56) |
             GT_SWAPBITPAIRS(kmer,56,4,52) |
             GT_SWAPBITPAIRS(kmer,54,6,48) |
             GT_SWAPBITPAIRS(kmer,52,8,44) |
             GT_SWAPBITPAIRS(kmer,50,10,40) |
             GT_SWAPBITPAIRS(kmer,48,12,36) |
             GT_SWAPBITPAIRS(kmer,46,14,32) |
             GT_SWAPBITPAIRS(kmer,44,16,28) |
             GT_SWAPBITPAIRS(kmer,42,18,24) |
             GT_SWAPBITPAIRS(kmer,40,20,20) |
             GT_SWAPBITPAIRS(kmer,38,22,16) |
             GT_SWAPBITPAIRS(kmer,36,24,12) |
             GT_SWAPBITPAIRS(kmer,34,26,8) |
             GT_SWAPBITPAIRS(kmer,32,28,4) |
             (kmer & (3U << 30));
    case 32:
      return GT_SWAPBITPAIRS(kmer,62,0,62) |
             GT_SWAPBITPAIRS(kmer,60,2,58) |
             GT_SWAPBITPAIRS(kmer,58,4,54) |
             GT_SWAPBITPAIRS(kmer,56,6,50) |
             GT_SWAPBITPAIRS(kmer,54,8,46) |
             GT_SWAPBITPAIRS(kmer,52,10,42) |
             GT_SWAPBITPAIRS(kmer,50,12,38) |
             GT_SWAPBITPAIRS(kmer,48,14,34) |
             GT_SWAPBITPAIRS(kmer,46,16,30) |
             GT_SWAPBITPAIRS(kmer,44,18,26) |
             GT_SWAPBITPAIRS(kmer,42,20,22) |
             GT_SWAPBITPAIRS(kmer,40,22,18) |
             GT_SWAPBITPAIRS(kmer,38,24,14) |
             GT_SWAPBITPAIRS(kmer,36,26,10) |
             GT_SWAPBITPAIRS(kmer,34,28,6) |
             GT_SWAPBITPAIRS(kmer,32,30,2);
#endif
    default: fprintf(stderr,"illegal kmersize=%u\n",kmersize);
             exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

static GtCodetype gt_kmercode_complement(GtCodetype kmer,GtCodetype maskright)
{
  return kmer ^ maskright;
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

#define UPDATEKMER(KMER,CC)\
        KMER <<= 2;\
        KMER |= CC;\
        KMER &= maskright

#define GT_ADJUSTREVERSEPOS(POS) (rightbound - (POS))

static GtCodetype getencseqkmers_nospecialtwobitencoding(
                                    const GtTwobitencoding *twobitencoding,
                                    unsigned long totallength,
                                    GtCodetype maskright,
                                    GtReadmode readmode,
                                    unsigned int kmersize,
                                    unsigned int upperkmersize,
                                    void(*processkmercode)(void *,
                                                           bool,
                                                           unsigned long,
                                                           GtCodetype),
                                    void *processkmercodeinfo,
                                    bool onlyfirst,
                                    unsigned long startpos,
                                    unsigned long endpos)
{
  unsigned long pos, unitindex, rightbound = totallength - kmersize;
  unsigned int shiftright;
  GtCodetype code;
  GtUchar cc;
  GtTwobitencoding currentencoding;

  gt_assert(kmersize > 1U);
  if (GT_ISDIRREVERSE(readmode))
  {
    unsigned long startpos2;

    gt_assert(endpos >= (unsigned long) upperkmersize);
    pos = endpos - (unsigned long) kmersize;
    unitindex = (pos > 0) ? GT_DIVBYUNITSIN2BITENC(pos-1) : 0;
    code = gt_kmercode_reverse(gt_kmercode_at_position(twobitencoding,pos,
                                                       kmersize),
                               kmersize);
    processkmercode(processkmercodeinfo,
                    true,
                    GT_ADJUSTREVERSEPOS(pos),
                    GT_ISDIRCOMPLEMENT(readmode)
                      ? gt_kmercode_complement(code,maskright)
                      : code);
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
      UPDATEKMER(code,cc);
      processkmercode(processkmercodeinfo,false,GT_ADJUSTREVERSEPOS(pos),
                       (readmode == GT_READMODE_REVCOMPL)
                          ? gt_kmercode_complement(code,maskright)
                          : code);
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
    unsigned long maxunitindex = gt_unitsoftwobitencoding(totallength) - 1;

    pos = startpos;
    unitindex = GT_DIVBYUNITSIN2BITENC(startpos+kmersize);
    code = gt_kmercode_at_position(twobitencoding,pos,kmersize);
    processkmercode(processkmercodeinfo,true,pos,
                    GT_ISDIRCOMPLEMENT(readmode)
                      ? gt_kmercode_complement(code,maskright)
                      : code);
    if (onlyfirst)
    {
      return code;
    }
    currentencoding = twobitencoding[unitindex];
    shiftright = (unsigned int)
                 GT_MULT2(GT_UNITSIN2BITENC - 1 -
                          GT_MODBYUNITSIN2BITENC(startpos+kmersize));
    while (pos < endpos - (unsigned long) upperkmersize)
    {
      pos++;
      cc = (GtUchar) (currentencoding >> shiftright) & 3;
      UPDATEKMER(code,cc);
      processkmercode(processkmercodeinfo,false,pos,
                       (readmode == GT_READMODE_COMPL)
                         ? gt_kmercode_complement(code,maskright)
                         : code);
      if (shiftright > 0)
      {
        shiftright -= 2;
      } else
      {
        gt_assert(unitindex < maxunitindex-1 ||
                  pos == endpos - (unsigned long) upperkmersize);
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

static void getencseqkmers_rangetwobitencoding(
                                      const GtTwobitencoding *twobitencoding,
                                      unsigned long totallength,
                                      unsigned long realtotallength,
                                      bool mirrored,
                                      GtCodetype maskright,
                                      GtReadmode readmode,
                                      unsigned int kmersize,
                                      unsigned int upperkmersize,
                                      bool onlyfirst,
                                      void(*processkmercode)(void *,
                                                             bool,
                                                             unsigned long,
                                                             GtCodetype),
                                      void *processkmercodeinfo,
                                      void(*processkmerspecial)(void *,
                                                                unsigned int,
                                                                unsigned int,
                                                                unsigned long),
                                      void *processkmerspecialinfo,
                                      unsigned long startpos,
                                      unsigned long endpos)
{
  GtCodetype lastcode, newcode;

  if (mirrored && startpos >= realtotallength) {
    gt_readmode_invert(readmode);
    startpos = GT_REVERSEPOS(realtotallength, startpos - realtotallength - 2);
    if (endpos == totallength)
      endpos = 0;
    else
      endpos = GT_REVERSEPOS(realtotallength, endpos - realtotallength - 2);
    if (startpos > endpos) {
      unsigned long tmp = startpos;
      startpos = endpos;
      endpos = tmp;
    }
    gt_assert(startpos <= endpos);
    gt_assert(endpos <= realtotallength);
  }
  if (endpos - startpos >= (unsigned long) upperkmersize)
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

      gt_assert((unsigned long) kmersize > endpos - startpos);
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

void getencseqkmers_twobitencoding(const GtEncseq *encseq,
                                   GtReadmode readmode,
                                   unsigned int kmersize,
                                   unsigned int upperkmersize,
                                   bool onlyfirst,
                                   void(*processkmercode)(void *,
                                                          bool,
                                                          unsigned long,
                                                          GtCodetype),
                                   void *processkmercodeinfo,
                                   void(*processkmerspecial)(void *,
                                                             unsigned int,
                                                             unsigned int,
                                                             unsigned long),
                                   void *processkmerspecialinfo)
{
  unsigned long laststart = 0, lastend,
                totallength,
                realtotallength;
  const GtTwobitencoding *twobitencoding
    = gt_encseq_twobitencoding_export(encseq);
  const GtCodetype maskright = GT_MASKRIGHT(kmersize);
  bool mirrored = gt_encseq_is_mirrored(encseq);

  lastend = totallength = realtotallength = gt_encseq_total_length(encseq);
  if (mirrored) {
    gt_assert((totallength & 1) == 1UL);
    realtotallength = ((realtotallength - 1) / 2);
  }
  if (gt_encseq_has_specialranges(encseq))
  {
    GtSpecialrangeiterator *sri;
    GtRange range;

    if (GT_ISDIRREVERSE(readmode))
    {
      sri = gt_specialrangeiterator_new(encseq,false);
      while (gt_specialrangeiterator_next(sri,&range))
      {
        gt_assert(range.end <= lastend);
        getencseqkmers_rangetwobitencoding(twobitencoding,
                                           totallength,
                                           realtotallength,
                                           mirrored,
                                           maskright,
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
    } else
    {
      sri = gt_specialrangeiterator_new(encseq,true);
      while (gt_specialrangeiterator_next(sri,&range))
      {
        gt_assert(range.start >= laststart);
        getencseqkmers_rangetwobitencoding(twobitencoding,
                                           totallength,
                                           realtotallength,
                                           mirrored,
                                           maskright,
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
    }
    gt_assert(totallength >= laststart);
    gt_specialrangeiterator_delete(sri);
  }
  getencseqkmers_rangetwobitencoding(twobitencoding,
                                     totallength,
                                     realtotallength,
                                     mirrored,
                                     maskright,
                                     readmode,
                                     kmersize,
                                     upperkmersize,
                                     onlyfirst,
                                     processkmercode,
                                     processkmercodeinfo,
                                     processkmerspecial,
                                     processkmerspecialinfo,
                                     GT_ISDIRREVERSE(readmode) ? 0
                                                               : laststart,
                                     GT_ISDIRREVERSE(readmode) ? lastend
                                                               : totallength);
}

static void gt_updateleftborderforkmer(Sfxiterator *sfi,
                                       GT_UNUSED bool firstinrange,
                                       GT_UNUSED unsigned long position,
                                       GtCodetype code)
{
  gt_assert(sfi->sfxstrategy.spmopt_minlength == 0);
  gt_bcktab_leftborder_addcode(sfi->leftborder,code);
}

static void gt_updateleftborderforspecialkmer(Sfxiterator *sfi,
                                              unsigned int maxprefixindex,
                                              unsigned int code,
                                              unsigned long position)
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

static void gt_spmopt_updateleftborderforkmer(Sfxiterator *sfi,
                                              bool firstinrange,
                                              GT_UNUSED unsigned long position,
                                              GtCodetype scancode)
{
  gt_assert(sfi->sfxstrategy.spmopt_minlength > 0);
  if (firstinrange || gt_checksuffixprefixbuckets(sfi,scancode))
  {
    gt_bcktab_leftborder_addcode(sfi->leftborder,
                                 GT_SCANCODE_TO_BCKCODE(scancode));
  }
}

static void gt_htinsertremainingcodes(GtUint64hashtable *table,
                                      bool firstinrange,
                                      GT_UNUSED unsigned long pos,
                                      GtCodetype code)
{
  if (!firstinrange)
  {
    (void) gt_uint64hashtable_search(table,(uint64_t) code,false);
  }
}

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

#define PROCESSKMERPREFIX(FUN) spmopt_updateleftborder_##FUN
#define PROCESSKMERTYPE        Sfxiterator
#define PROCESSKMERSPECIALTYPE GT_UNUSED Sfxiterator
#define PROCESSKMERCODE        gt_spmopt_updateleftborderforkmer

#include "sfx-mapped4.gen"

#undef PROCESSKMERPREFIX
#undef PROCESSKMERTYPE
#undef PROCESSKMERSPECIALTYPE
#undef PROCESSKMERCODE

#define PROCESSKMERPREFIX(FUN)          insertsuffix_##FUN
#define PROCESSKMERTYPE                 void
#define PROCESSKMERSPECIALTYPE          GT_UNUSED Sfxiterator
#define PROCESSKMERCODE                 gt_insertkmerwithoutspecial1

#include "sfx-mapped4.gen"

#undef PROCESSKMERPREFIX
#undef PROCESSKMERTYPE
#undef PROCESSKMERSPECIALTYPE
#undef PROCESSKMERCODE

#define PROCESSKMERPREFIX(FUN)          htinsertsuffixremainingcodes_##FUN
#define PROCESSKMERTYPE                 GtUint64hashtable
#define PROCESSKMERSPECIALTYPE          GT_UNUSED void
#define PROCESSKMERCODE                 gt_htinsertremainingcodes

#define GT_MAPPED4_GLOBAL
#include "sfx-mapped4.gen"
#undef GT_MAPPED4_GLOBAL

/*
#define SHOWCURRENTSPACE\
        printf("spacepeak at line %d: %.2f\n",__LINE__,\
          GT_MEGABYTES(gt_ma_get_space_current() + gt_fa_get_space_current()))
*/

#define SHOWCURRENTSPACE /* Nothing */

static void gt_sfimarkprefixsuffixbuckets(void *processinfo,
                                          GT_UNUSED bool firstinrange,
                                          GT_UNUSED unsigned long pos,
                                          GtCodetype scancode)
{
  Sfxiterator *sfi = (Sfxiterator *) processinfo;
  GtCodetype checkcode = GT_SCANCODE_TO_PREFIXCODE(scancode);

  gt_assert(firstinrange);
  if (!GT_ISIBITSET(sfi->markprefixbuckets,checkcode))
  {
    GT_SETIBIT(sfi->markprefixbuckets,checkcode);
  }
#ifdef _LP64
  checkcode = GT_SCANCODE_TO_SUFFIXCODE(scancode);
  if (!GT_ISIBITSET(sfi->marksuffixbuckets,checkcode))
  {
    GT_SETIBIT(sfi->marksuffixbuckets,checkcode);
  }
#endif
}

static size_t gt_sizeforbittable(unsigned int numofchars,
                                 unsigned int prefixlength)
{
  unsigned long numofcodes;

  numofcodes = gt_power_for_small_exponents(numofchars,prefixlength);
  return sizeof (GtBitsequence) * GT_NUMOFINTSFORBITS(numofcodes);
}

static void gt_determineaddionalsuffixprefixchars(
                                   unsigned int *additionalprefixchars,
                                   unsigned int *additionalsuffixchars,
                                   unsigned int numofchars,
                                   unsigned int prefixlength,
                                   size_t estimatedspace,
                                   unsigned long maximumspace)
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
#ifdef _LP64
  {
    unsigned int suffixchars;
    size_t sizeofsuffixmarks;
    sizeofprefixmarks = gt_sizeforbittable(numofchars,prefixlength+prefixchars);
    for (suffixchars = 1U;
         prefixlength+prefixchars+prefixlength+suffixchars <= GT_UNITSIN2BITENC;
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

static unsigned long gt_bcktab_code_to_prefix_index(unsigned long code,
                                             unsigned int additionalprefixchars)
{
  if (GT_MULT2(additionalprefixchars) > (unsigned int) GT_LOGWORDSIZE)
  {
    return (unsigned long) (code << (GT_MULT2(additionalprefixchars) -
                                 GT_LOGWORDSIZE));
  }
  return (unsigned long) (code >> (GT_LOGWORDSIZE -
                                GT_MULT2(additionalprefixchars)));
}

Sfxiterator *gt_Sfxiterator_new_withadditionalvalues(
                                const GtEncseq *encseq,
                                GtReadmode readmode,
                                unsigned int prefixlength,
                                unsigned int numofparts,
                                unsigned long maximumspace,
                                void *voidoutlcpinfo,
                                FILE *outfpbcktab,
                                const Sfxstrategy *sfxstrategy,
                                GtTimer *sfxprogress,
                                bool withprogressbar,
                                GtLogger *logger,
                                GtError *err)
{
  Sfxiterator *sfi = NULL;
  unsigned long realspecialranges, specialcharacters, numofsuffixestosort = 0;
  bool haserr = false;
#ifdef _LP64
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
        prefixlength > (unsigned int) MAXPREFIXLENGTH)
    {
      gt_error_set(err,"argument for option -pl must be in the range [1,%u]",
                    MAXPREFIXLENGTH);
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
      gt_logger_log(logger,"sizeof (spaceCodeatposition)=%lu bytes",
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
    sfi->totallength = gt_encseq_total_length(encseq);
    gt_logger_log(logger,"totallength=%lu",sfi->totallength);
    sfi->specialcharacters = specialcharacters;
    sfi->outlcpinfo = (Outlcpinfo *) voidoutlcpinfo;
    sfi->outlcpinfoforsample = NULL;
    sfi->sri = NULL;
    sfi->part = 0;
    sfi->exhausted = false;
    sfi->bucketiterstep = 0;
    sfi->logger = logger;
    sfi->sfxprogress = sfxprogress;
    if (sfi->sfxstrategy.differencecover > 0 &&
        specialcharacters < sfi->totallength)
    {
      if (sfi->outlcpinfo != NULL)
      {
        sfi->outlcpinfoforsample
          = gt_Outlcpinfo_new(NULL,sfi->numofchars,0,err);
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
          estimatedspace += gt_differencecover_requiredspace(sfi->dcov);
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
    if (prefixlength > 1U && gt_has_twobitencoding(sfi->encseq) &&
        sfi->sfxstrategy.spmopt_minlength > 0)
    {
      unsigned int suffixchars = 0, additionalsuffixchars = 2U;
      size_t sizeofprefixmarks, intsforbits;
#ifdef _LP64
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
#ifdef _LP64
      if (sfi->prefixlength + sfi->spmopt_additionalprefixchars +
          sfi->prefixlength + additionalsuffixchars >
          (unsigned int) GT_UNITSIN2BITENC)
      {
        suffixchars = GT_UNITSIN2BITENC -
                      (sfi->prefixlength + sfi->spmopt_additionalprefixchars);
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
                                "table of %lu bytes",
                                sfi->prefixlength +
                                sfi->spmopt_additionalprefixchars,
                                (unsigned long) sizeofprefixmarks);
      sfi->mappedmarkprefixbuckets
        = gt_Sfxmappedrange_new_basic(sfi->spmopt_numofallprefixcodes,
                                      GtSfxGtBitsequence,
                                      gt_bcktab_code_to_prefix_index,
                                      sfi->spmopt_additionalprefixchars);
      sfi->spmopt_numofallsuffixcodes
        = gt_power_for_small_exponents(sfi->numofchars,suffixchars);
#ifdef _LP64
      GT_INITBITTAB(sfi->marksuffixbuckets,
                    sfi->spmopt_numofallsuffixcodes);
      sizeofsuffixmarks
        = sizeof (*sfi->marksuffixbuckets) *
                  GT_NUMOFINTSFORBITS(sfi->spmopt_numofallsuffixcodes);
      estimatedspace += sizeofsuffixmarks;
      gt_logger_log(sfi->logger,"for all sequences, keep track of "
                                "%u-mers starting at position %u "
                                "using a table of %lu bytes",
                                suffixchars,
                                sfi->prefixlength +
                                sfi->spmopt_additionalprefixchars,
                                (unsigned long) sizeofsuffixmarks);
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
    unsigned long largestbucketsize,
                  saved_bucketswithoutwholeleaf;
    gt_assert(sfi != NULL);
    sfi->storespecials = true;
    if (sfxprogress != NULL)
    {
      gt_timer_show_progress(sfxprogress, "counting prefix distribution",
                             stdout);
    }
    if (prefixlength == 1U)
    {
      unsigned int charidx;

      for (charidx=0; charidx<sfi->numofchars; charidx++)
      {
        unsigned int updateindex = GT_ISDIRCOMPLEMENT(readmode) ?
                                        GT_COMPLEMENTBASE(charidx) :
                                        charidx;
        gt_bcktab_leftborder_assign(sfi->leftborder,(GtCodetype) updateindex,
                                    gt_encseq_charcount(encseq,
                                                        (GtUchar) charidx));
      }
    } else
    {
      if (gt_has_twobitencoding(encseq) &&
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
        gt_assert(realspecialranges+1 >= sfi->nextfreeCodeatposition);
        reversespecialcodes(sfi->spaceCodeatposition,
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
    gt_logger_log(sfi->logger, "largest bucket size=%lu",largestbucketsize);
    if (sfi->sfxstrategy.spmopt_minlength > 0)
    {
      gt_logger_log(sfi->logger, "relevant suffixes=%.2f%%",100.0 *
                                        (double) numofsuffixestosort/
                                        (sfi->totallength+1));
      gt_logger_log(sfi->logger,"saved_bucketswithoutwholeleaf=%lu",
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
    if (maximumspace > 0)
    {
      int retval;

      gt_assert(numofparts == 1U);
      retval = computepartsfittingmaximumspace(
                                       estimatedspace,
                                       maximumspace,
                                       sfi->bcktab,
                                       sfi->mappedmarkprefixbuckets,
                                       sfi->totallength,
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
    sfi->suftabparts = gt_newsuftabparts(numofparts,
                                         sfi->bcktab,
                                         sfi->mappedmarkprefixbuckets,
                                         numofsuffixestosort,
                                         specialcharacters + 1,
                                         logger);
    gt_assert(sfi->suftabparts != NULL);
    if (stpgetnumofparts(sfi->suftabparts) > 1U)
    {
      if (gt_bcktab_storetmp(sfi->bcktab, sfi->logger, err) != 0)
      {
        haserr = true;
      }
    }
  }
  SHOWACTUALSPACE;
  if (!haserr)
  {
    gt_assert(sfi != NULL && sfi->suftabparts != NULL);
    if (stpgetnumofparts(sfi->suftabparts) > 1U &&
        sfi->sfxstrategy.spmopt_minlength > 0)
    {
      gt_assert(sfi->markprefixbuckets != NULL);
      gt_assert(sfi->mappedmarkprefixbuckets != NULL);
      if (gt_Sfxmappedrange_enhance(sfi->mappedmarkprefixbuckets,
                                    (void **) &sfi->markprefixbuckets,
                                    false,
                                    "markprefixbuckets",
                                    sfi->logger,
                                    err) != 0)
      {
        sfi->markprefixbuckets = NULL;
        haserr = true;
      }
      gt_assert(sfi->markprefixbuckets == NULL);
    } else
    {
      gt_Sfxmappedrange_delete(sfi->mappedmarkprefixbuckets,logger);
      sfi->mappedmarkprefixbuckets = NULL;
    }
  }
  SHOWACTUALSPACE;
  if (!haserr)
  {
    gt_assert(sfi != NULL);
    sfi->suffixsortspace
      = gt_suffixsortspace_new(stpgetlargestsuftabwidth(sfi->suftabparts),
                               sfi->totallength,
                               sfi->sfxstrategy.suftabuint,
                               logger);
    gt_assert(sfi->suffixsortspace);
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
    sfi->fusp.allocatedSuffixptr = stpgetlargestsuftabwidth(sfi->suftabparts);
    sfi->overhang.start = sfi->overhang.end = 0;
  }
  SHOWACTUALSPACE;
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
                                unsigned long maximumspace,
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

static void gt_sfxiterator_preparethispart(Sfxiterator *sfi)
{
  unsigned long partwidth;
  GtBucketspec2 *bucketspec2 = NULL;

  if (sfi->part == 0 && sfi->withprogressbar)
  {
    gt_assert(sfi->bcktab != NULL);
    gt_progressbar_start(&sfi->bucketiterstep,
                         (unsigned long long)
                         gt_bcktab_numofallcodes(sfi->bcktab));
  }
  sfi->currentmincode = stpgetcurrentmincode(sfi->part,sfi->suftabparts);
  sfi->currentmaxcode = stpgetcurrentmaxcode(sfi->part,sfi->suftabparts);
  sfi->widthofpart = stpgetcurrentwidthofpart(sfi->part,sfi->suftabparts);
  if (sfi->sfxprogress != NULL)
  {
    gt_timer_show_progress(sfi->sfxprogress, "inserting suffixes into buckets",
                           stdout);
  }
  if (stpgetnumofparts(sfi->suftabparts) > 1U)
  {
    gt_logger_log(sfi->logger,"compute part %u: "
                              "%lu suffixes,%lu buckets from "
                              "%lu..%lu",
                  sfi->part,
                  stpgetcurrentwidthofpart(sfi->part,sfi->suftabparts),
                  sfi->currentmaxcode - sfi->currentmincode + 1,
                  sfi->currentmincode,
                  sfi->currentmaxcode);
    gt_bcktab_assignboundsforpart(sfi->bcktab,
                                  sfi->part,
                                  sfi->currentmincode,
                                  sfi->currentmaxcode,
                                  sfi->logger);
    if (sfi->mappedmarkprefixbuckets != NULL)
    {
      sfi->markprefixbuckets
        = (GtBitsequence *)
          gt_Sfxmappedrange_map(sfi->mappedmarkprefixbuckets,
                                sfi->part,
                                sfi->currentmincode,
                                sfi->currentmaxcode,
                                sfi->logger);
    }
  }
  SHOWACTUALSPACE;
  gt_suffixsortspace_partoffset_set(sfi->suffixsortspace,
                                    stpgetcurrentsuftaboffset(sfi->part,
                                                             sfi->suftabparts));
  if (sfi->sfxstrategy.spmopt_minlength == 0)
  {
    if (sfi->sfxstrategy.storespecialcodes)
    {
      sfx_derivespecialcodesfromtable(sfi,
                    stpgetnumofparts(sfi->suftabparts) == 1U ? true : false);
    } else
    {
      sfx_derivespecialcodesonthefly(sfi);
    }
  }
  SHOWACTUALSPACE;
  sfi->exportptr = gt_suffixsortspace_exportptr(0,sfi->suffixsortspace);
  if (sfi->prefixlength > 1U
      && gt_has_twobitencoding(sfi->encseq)
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
  partwidth = stpgetcurrentsumofwdith(sfi->part,sfi->suftabparts);

  if (stpgetnumofparts(sfi->suftabparts) == 1U && sfi->outlcpinfo == NULL &&
      sfi->prefixlength >= 2U && sfi->sfxstrategy.spmopt_minlength == 0)
  {
    bucketspec2 = gt_copysort_new(sfi->bcktab,sfi->encseq,sfi->readmode,
                                  partwidth,sfi->numofchars);
  }
  SHOWACTUALSPACE;
  if (sfi->sfxstrategy.differencecover > 0)
  {
    gt_differencecoversetsuffixsortspace(sfi->dcov,sfi->suffixsortspace);
  }
  if (sfi->part == 0)
  {
    gt_logger_log(sfi->logger,"used workspace for sorting: %.2f MB",
                    GT_MEGABYTES(gt_size_of_sort_workspace (&sfi->sfxstrategy,
                                                            sfi->encseq,
                                                            sfi->prefixlength,
                                                            sfi->readmode)));
  }
  if (!sfi->sfxstrategy.onlybucketinsertion)
  {
    gt_sortallbuckets(sfi->suffixsortspace,
                      partwidth,
                      bucketspec2,
                      sfi->encseq,
                      sfi->readmode,
                      sfi->currentmincode,
                      sfi->currentmaxcode,
                      sfi->bcktab,
                      sfi->numofchars,
                      sfi->prefixlength,
                      sfi->sfxstrategy.differencecover == 0
                        ? sfi->outlcpinfo : NULL,
                      sfi->sfxstrategy.differencecover,
                      &sfi->sfxstrategy,
                      sfi->sfxstrategy.differencecover == 0
                        ? NULL : gt_differencecover_sortunsortedbucket,
                      sfi->sfxstrategy.differencecover == 0
                        ? NULL : (void *) sfi->dcov,
                      &sfi->bucketiterstep,
                      sfi->logger);
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
      gt_suffixsortspace_setdirect(sfi->fusp.sssp,
                                   sfi->fusp.nextfreeSuffixptr,
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
                                   sfi->fusp.nextfreeSuffixptr,
                                   pos);
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
        insertfullspecialrange(sfi,sfi->overhang.start,sfi->overhang.end);
        sfi->overhang.start = sfi->overhang.end = 0;
        break;
      }
      /* overhang fits into the buffer and buffer is not full */
      insertfullspecialrange(sfi,sfi->overhang.start,sfi->overhang.end);
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
            insertfullspecialrange(sfi,range.start + rest, range.end);
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
  sfi->fusp.nextfreeSuffixptr = 0;
  if (sfi->sfxstrategy.spmopt_minlength == 0)
  {
    gt_suffixsortspace_partoffset_set(sfi->suffixsortspace,
                                      sfi->totallength-sfi->specialcharacters);
    fillspecialnextpage(sfi);
    gt_assert(sfi->fusp.nextfreeSuffixptr > 0);
  } else
  {
    sfi->exhausted = true;
  }
  *numberofsuffixes = sfi->fusp.nextfreeSuffixptr;
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
  if (stpgetnumofparts(sfi->suftabparts) <= 1U)
  {
    int ret = gt_bcktab_flush_to_file(fp,sfi->bcktab,err);
    gt_fa_fclose(fp);
    return ret;
  }
  return 0;
}

unsigned long gt_Sfxiterator_longest(const Sfxiterator *sfi)
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
