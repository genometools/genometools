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
  GtCodetype numofallcodes;
  GtLeftborder *leftborder; /* points to bcktab->leftborder */
  GtBitsequence *markwholeleafbuckets;
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
  bool storespecials;
  unsigned long nextfreeCodeatposition;
  Codeatposition *spaceCodeatposition;
  GtStr *bcktmpfilename;
  unsigned int kmerfastmaskright;
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

  if (kmercode->definedspecialposition)
  {
    if (sfi->sfxstrategy.spmopt == 0)
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
    if (sfi->markwholeleafbuckets != NULL &&
        (position == 0 || gt_encseq_position_is_separator(sfi->encseq,
                                                          position - 1,
                                                          sfi->readmode)))
    {
      GT_SETIBIT(sfi->markwholeleafbuckets,kmercode->code);
    }
  }
#ifdef SKDEBUG
  previouscode = kmercode->code;
  previouskmercodedefined = kmercode->defined;
  previousstorespecials = sfi->storespecials;
  previousspecialpos = kmercode->specialposition;
#endif
}

static void gt_insertkmerwithoutspecial1(void *processinfo,
                                         GT_UNUSED bool firstinrange,
                                         unsigned long position,
                                         GtCodetype code)
{
  Sfxiterator *sfi = (Sfxiterator *) processinfo;

  if (code >= sfi->currentmincode && code <= sfi->currentmaxcode &&
      (sfi->markwholeleafbuckets == NULL ||
       GT_ISIBITSET(sfi->markwholeleafbuckets,code)))
  {
    unsigned long stidx;

    stidx = gt_bcktab_leftborder_insertionindex(sfi->leftborder,code);
    gt_suffixsortspace_setdirectwithoffset(sfi->suffixsortspace,stidx,
                                           position);
    /* from right to left */
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
          gt_bcktab_updatespecials(sfi->bcktab,code,sfi->numofchars,
                                   prefixindex);
          stidx = gt_bcktab_leftborder_insertionindex(sfi->leftborder,code);
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
            gt_bcktab_updatespecials(sfi->bcktab,code,sfi->numofchars,
                                     prefixindex);
            gt_assert(code > 0);
            stidx = gt_bcktab_leftborder_insertionindex(sfi->leftborder,code);
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
  if (sfi->bcktab != NULL)
  {
    gt_bcktab_addfinalspecials(sfi->bcktab,sfi->numofchars,
                               sfi->specialcharacters);
  }
  if (sfi->sri != NULL)
  {
    gt_specialrangeiterator_delete(sfi->sri);
  }
  gt_free(sfi->spaceCodeatposition);
  sfi->spaceCodeatposition = NULL;
  gt_suffixsortspace_delete(sfi->suffixsortspace,
                            sfi->sfxstrategy.spmopt == 0 ? true : false);
  gt_freesuftabparts(sfi->suftabparts);
  gt_bcktab_delete(sfi->bcktab);
  gt_Outlcpinfo_delete(sfi->outlcpinfoforsample);
  gt_free(sfi->markwholeleafbuckets);
  gt_differencecover_delete(sfi->dcov);
  if (sfi->bcktmpfilename != NULL)
  {
    gt_logger_log(sfi->logger,"remove \"%s\"",gt_str_get(sfi->bcktmpfilename));
    if (unlink(gt_str_get(sfi->bcktmpfilename)) != 0)
    {
      if (err != NULL)
      {
        gt_error_set(err,"Cannot unlink file \"%s\": %s",
                        gt_str_get(sfi->bcktmpfilename),
                        strerror(errno));
        haserr = true;
      } else
      {
        fprintf(stderr,"Cannot unlink file \"%s\": %s",
                        gt_str_get(sfi->bcktmpfilename),
                        strerror(errno));
        exit(EXIT_FAILURE);
      }
    }
  }
  gt_str_delete(sfi->bcktmpfilename);
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
                                           unsigned long totallength,
                                           unsigned long specialcharacters,
                                           unsigned long numofsuffixestosort,
                                           bool suftabcompressedbytes,
                                           GtError *err)
{
  unsigned int parts;
  Suftabparts *suftabparts;

  if (estimatedspace >= (size_t) maximumspace)
  {
    gt_error_set(err,"already used %.2f MB of memory, cannot compute "
                     "enhanced suffix array in at most %lu bytes",
                     GT_MEGABYTES(estimatedspace), maximumspace);
    return -1;
  }
  for (parts = 1U; parts <= 500U; parts++)
  {
    size_t suftabsize;

    suftabparts = gt_newsuftabparts(parts,
                                    bcktab,
                                    numofsuffixestosort,
                                    specialcharacters + 1,
                                    NULL);
    gt_assert(suftabparts != NULL);
    suftabsize = gt_suffixsortspace_requiredspace(
                                              stpgetlargestwidth(suftabparts),
                                              totallength,
                                              suftabcompressedbytes);
    if ((unsigned long) (suftabsize + estimatedspace) <= maximumspace)
    {
      gt_freesuftabparts(suftabparts);
      return (int) parts;
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
  unsigned int unitoffset = (unsigned int) GT_MODBYUNITSIN2BITENC(pos);
  unsigned long unitindex = GT_DIVBYUNITSIN2BITENC(pos);
  const GtCodetype maskright = (GtCodetype) (1 << GT_MULT2(kmersize))-1;

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
  const GtCodetype maskright = (GtCodetype) (1 << GT_MULT2(kmersize))-1;
  return (GtCodetype) (twobitencoding[0] >>
                       GT_MULT2(GT_UNITSIN2BITENC - kmersize)) & maskright;
}

#define GT_SWAPBITPAIRS(L1,L2,D) ((kmer & (3UL << L1)) >> D) |\
                                 ((kmer & (3UL << L2)) << D)

GtCodetype gt_kmercode_reverse(GtCodetype kmer,unsigned int kmersize)
{
  switch (kmersize)
  {
    case 2:
      return GT_SWAPBITPAIRS(2,0,2);
    case 3:
      return GT_SWAPBITPAIRS(4,0,4) |
             (kmer & (3U << 2));
    case 4:
      return GT_SWAPBITPAIRS(6,0,6) |
             GT_SWAPBITPAIRS(4,2,2);
    case 5:
      return GT_SWAPBITPAIRS(8,0,8) |
             GT_SWAPBITPAIRS(6,2,4) |
             (kmer & (3U << 4));
    case 6:
      return GT_SWAPBITPAIRS(10,0,10) |
             GT_SWAPBITPAIRS(8,2,6) |
             GT_SWAPBITPAIRS(6,4,2);
    case 7:
      return GT_SWAPBITPAIRS(12,0,12) |
             GT_SWAPBITPAIRS(10,2,8) |
             GT_SWAPBITPAIRS(8,4,4) |
             (kmer & (3U << 6));
    case 8:
      return GT_SWAPBITPAIRS(14,0,14) |
             GT_SWAPBITPAIRS(12,2,10) |
             GT_SWAPBITPAIRS(10,4,6) |
             GT_SWAPBITPAIRS(8,6,2);
    case 9:
      return GT_SWAPBITPAIRS(16,0,16) |
             GT_SWAPBITPAIRS(14,2,12) |
             GT_SWAPBITPAIRS(12,4,8) |
             GT_SWAPBITPAIRS(10,6,4) |
             (kmer & (3U << 8));
    case 10:
      return GT_SWAPBITPAIRS(18,0,18) |
             GT_SWAPBITPAIRS(16,2,14) |
             GT_SWAPBITPAIRS(14,4,10) |
             GT_SWAPBITPAIRS(12,6,6) |
             GT_SWAPBITPAIRS(10,8,2);
    case 11:
      return GT_SWAPBITPAIRS(20,0,20) |
             GT_SWAPBITPAIRS(18,2,16) |
             GT_SWAPBITPAIRS(16,4,12) |
             GT_SWAPBITPAIRS(14,6,8) |
             GT_SWAPBITPAIRS(12,8,4) |
             (kmer & (3U << 10));
    case 12:
      return GT_SWAPBITPAIRS(22,0,22) |
             GT_SWAPBITPAIRS(20,2,18) |
             GT_SWAPBITPAIRS(18,4,14) |
             GT_SWAPBITPAIRS(16,6,10) |
             GT_SWAPBITPAIRS(14,8,6) |
             GT_SWAPBITPAIRS(12,10,2);
    case 13:
      return GT_SWAPBITPAIRS(24,0,24) |
             GT_SWAPBITPAIRS(22,2,20) |
             GT_SWAPBITPAIRS(20,4,16) |
             GT_SWAPBITPAIRS(18,6,12) |
             GT_SWAPBITPAIRS(16,8,8) |
             GT_SWAPBITPAIRS(14,10,4) |
             (kmer & (3U << 12));
    case 14:
      return GT_SWAPBITPAIRS(26,0,26) |
             GT_SWAPBITPAIRS(24,2,22) |
             GT_SWAPBITPAIRS(22,4,18) |
             GT_SWAPBITPAIRS(20,6,14) |
             GT_SWAPBITPAIRS(18,8,10) |
             GT_SWAPBITPAIRS(16,10,6) |
             GT_SWAPBITPAIRS(14,12,2);
#ifdef _LP64
    case 15:
      return GT_SWAPBITPAIRS(28,0,28) |
             GT_SWAPBITPAIRS(26,2,24) |
             GT_SWAPBITPAIRS(24,4,20) |
             GT_SWAPBITPAIRS(22,6,16) |
             GT_SWAPBITPAIRS(20,8,12) |
             GT_SWAPBITPAIRS(18,10,8) |
             GT_SWAPBITPAIRS(16,12,4) |
             (kmer & (3U << 14));
    case 16:
      return GT_SWAPBITPAIRS(30,0,30) |
             GT_SWAPBITPAIRS(28,2,26) |
             GT_SWAPBITPAIRS(26,4,22) |
             GT_SWAPBITPAIRS(24,6,18) |
             GT_SWAPBITPAIRS(22,8,14) |
             GT_SWAPBITPAIRS(20,10,10) |
             GT_SWAPBITPAIRS(18,12,6) |
             GT_SWAPBITPAIRS(16,14,2);
    case 17:
      return GT_SWAPBITPAIRS(32,0,32) |
             GT_SWAPBITPAIRS(30,2,28) |
             GT_SWAPBITPAIRS(28,4,24) |
             GT_SWAPBITPAIRS(26,6,20) |
             GT_SWAPBITPAIRS(24,8,16) |
             GT_SWAPBITPAIRS(22,10,12) |
             GT_SWAPBITPAIRS(20,12,8) |
             GT_SWAPBITPAIRS(18,14,4) |
             (kmer & (3U << 16));
    case 18:
      return GT_SWAPBITPAIRS(34,0,34) |
             GT_SWAPBITPAIRS(32,2,30) |
             GT_SWAPBITPAIRS(30,4,26) |
             GT_SWAPBITPAIRS(28,6,22) |
             GT_SWAPBITPAIRS(26,8,18) |
             GT_SWAPBITPAIRS(24,10,14) |
             GT_SWAPBITPAIRS(22,12,10) |
             GT_SWAPBITPAIRS(20,14,6) |
             GT_SWAPBITPAIRS(18,16,2);
    case 19:
      return GT_SWAPBITPAIRS(36,0,36) |
             GT_SWAPBITPAIRS(34,2,32) |
             GT_SWAPBITPAIRS(32,4,28) |
             GT_SWAPBITPAIRS(30,6,24) |
             GT_SWAPBITPAIRS(28,8,20) |
             GT_SWAPBITPAIRS(26,10,16) |
             GT_SWAPBITPAIRS(24,12,12) |
             GT_SWAPBITPAIRS(22,14,8) |
             GT_SWAPBITPAIRS(20,16,4) |
             (kmer & (3U << 18));
    case 20:
      return GT_SWAPBITPAIRS(38,0,38) |
             GT_SWAPBITPAIRS(36,2,34) |
             GT_SWAPBITPAIRS(34,4,30) |
             GT_SWAPBITPAIRS(32,6,26) |
             GT_SWAPBITPAIRS(30,8,22) |
             GT_SWAPBITPAIRS(28,10,18) |
             GT_SWAPBITPAIRS(26,12,14) |
             GT_SWAPBITPAIRS(24,14,10) |
             GT_SWAPBITPAIRS(22,16,6) |
             GT_SWAPBITPAIRS(20,18,2);
#endif
    default: fprintf(stderr,"%s: illegal kmersize=%u\n",__func__,kmersize);
             exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

GtCodetype gt_kmercode_complement(GtCodetype kmer,GtCodetype maskright)
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

#define ADJUSTREVERSEPOS(POS) (rightbound - (POS))

static GtCodetype getencseqkmers_nospecialtwobitencoding(
                                    const GtTwobitencoding *twobitencoding,
                                    unsigned long totallength,
                                    GtCodetype maskright,
                                    GtReadmode readmode,
                                    unsigned int kmersize,
                                    void(*processkmercode)(void *,
                                                           unsigned long,
                                                           GtCodetype),
                                    void *processkmercodeinfo,
                                    unsigned long startpos,
                                    unsigned long endpos)
{
  unsigned long pos, unitindex;
  unsigned int shiftright;
  GtCodetype code;
  GtUchar cc;
  GtTwobitencoding currentencoding;
  unsigned long rightbound = totallength - kmersize;

  gt_assert(kmersize > 1U);
  if (GT_ISDIRREVERSE(readmode))
  {
    gt_assert(endpos >= (unsigned long) kmersize);
    pos = endpos - (unsigned long) kmersize;
    unitindex = (pos > 0) ? GT_DIVBYUNITSIN2BITENC(pos-1) : 0;
    code = gt_kmercode_reverse(gt_kmercode_at_position(twobitencoding,pos,
                                                       kmersize),
                               kmersize);
    processkmercode(processkmercodeinfo,
                    ADJUSTREVERSEPOS(pos),
                    GT_ISDIRCOMPLEMENT(readmode)
                      ? gt_kmercode_complement(code,maskright)
                      : code);
  } else
  {
    pos = startpos;
    unitindex = GT_DIVBYUNITSIN2BITENC(startpos+kmersize);
    code = gt_kmercode_at_position(twobitencoding,pos,kmersize);
    processkmercode(processkmercodeinfo,pos,
                    GT_ISDIRCOMPLEMENT(readmode)
                      ? gt_kmercode_complement(code,maskright)
                      : code);
  }
  currentencoding = twobitencoding[unitindex];
  if (GT_ISDIRREVERSE(readmode))
  {
    shiftright = (unsigned int)
                 GT_MULT2(GT_UNITSIN2BITENC - 1 -
                          GT_MODBYUNITSIN2BITENC(pos-1));
    while (pos > startpos)
    {
      pos--;
      cc = (GtUchar) (currentencoding >> shiftright) & 3;
      UPDATEKMER(code,cc);
      processkmercode(processkmercodeinfo,ADJUSTREVERSEPOS(pos),
                       (readmode == GT_READMODE_REVCOMPL)
                          ? gt_kmercode_complement(code,maskright)
                          : code);
      if (shiftright < (unsigned int) (GT_INTWORDSIZE-2))
      {
        shiftright += 2;
      } else
      {
        gt_assert(unitindex > 0 || pos == startpos);
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

    shiftright = (unsigned int)
                 GT_MULT2(GT_UNITSIN2BITENC - 1 -
                          GT_MODBYUNITSIN2BITENC(startpos+kmersize));
    while (pos < endpos - (unsigned long) kmersize)
    {
      pos++;
      cc = (GtUchar) (currentencoding >> shiftright) & 3;
      UPDATEKMER(code,cc);
      processkmercode(processkmercodeinfo,pos,
                       (readmode == GT_READMODE_COMPL)
                         ? gt_kmercode_complement(code,maskright)
                         : code);
      if (shiftright > 0)
      {
        shiftright -= 2;
      } else
      {
        gt_assert(unitindex < maxunitindex-1 ||
                  pos == endpos - (unsigned long) kmersize);
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
                                      void(*processkmercode)(void *,
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
  if (endpos - startpos >= (unsigned long) kmersize)
  {
    gt_assert(endpos > 0);
    lastcode = getencseqkmers_nospecialtwobitencoding(twobitencoding,
                                                      totallength,
                                                      maskright,
                                                      readmode,
                                                      kmersize,
                                                      processkmercode,
                                                      processkmercodeinfo,
                                                      startpos,
                                                      endpos);
    if (GT_ISDIRCOMPLEMENT(readmode))
    {
      lastcode = gt_kmercode_complement(lastcode,maskright);
    }
    newcode = ((lastcode << 2) | 3UL) & maskright;
    if (processkmerspecial != NULL)
    {
      processkmerspecial(processkmerspecialinfo,
                         kmersize-1,
                         (unsigned int) newcode,
                         GT_ISDIRREVERSE(readmode) ? (totallength - startpos)
                                                   : endpos);
    }
  } else
  {
    if (startpos < endpos)
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
      if (processkmerspecial != NULL)
      {
        processkmerspecial(processkmerspecialinfo,
                           (unsigned int) (endpos - startpos),
                           (unsigned int) newcode,
                           GT_ISDIRREVERSE(readmode) ? (totallength - startpos)
                                                     : endpos);
      }
    }
  }
}

void getencseqkmers_twobitencoding(const GtEncseq *encseq,
                                   GtReadmode readmode,
                                   unsigned int kmersize,
                                   void(*processkmercode)(void *,
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
  const GtCodetype maskright = (GtCodetype) (1 << GT_MULT2(kmersize))-1;
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
  gt_assert(sfi->sfxstrategy.spmopt == 0);
  gt_bcktab_leftborder_addcode(sfi->leftborder,code);
}

static void gt_updateleftborderforspecialkmer(Sfxiterator *sfi,
                                              unsigned int maxprefixindex,
                                              unsigned int code,
                                              unsigned long position)
{
  unsigned int idx;

  gt_assert(sfi->sfxstrategy.spmopt == 0);
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
                                              GtCodetype code)
{
  gt_assert(sfi->sfxstrategy.spmopt > 0);
  gt_bcktab_leftborder_addcode(sfi->leftborder,code);
  if (firstinrange)
  {
    GT_SETIBIT(sfi->markwholeleafbuckets,code);
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

/*
#define SHOWCURRENTSPACE\
        printf("spacepeak at line %d: %.2f\n",__LINE__,\
          GT_MEGABYTES(gt_ma_get_space_current() + gt_fa_get_space_current()))
*/

#define SHOWCURRENTSPACE /* Nothing */

Sfxiterator *gt_Sfxiterator_new(const GtEncseq *encseq,
                                GtReadmode readmode,
                                unsigned int prefixlength,
                                unsigned int numofparts,
                                unsigned long maximumspace,
                                void *voidoutlcpinfo,
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
      if (sfxstrategy->spmopt > 0 && prefixlength > sfxstrategy->spmopt)
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
        sfxstrategy->spmopt == 0)
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
    sfi->numofchars = gt_encseq_alphabetnumofchars(encseq);
    sfi->prefixlength = prefixlength;
    sfi->kmerfastmaskright = (1U << GT_MULT2(prefixlength))-1;
    sfi->dcov = NULL;
    sfi->markwholeleafbuckets = NULL;
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
  SHOWCURRENTSPACE;
  if (!haserr)
  {
    uint64_t sizeofbcktab;
    gt_assert(sfi != NULL);
    sfi->bcktab = gt_bcktab_new(sfi->numofchars,
                                prefixlength,
                                sfi->totallength+1,
                                sfi->sfxstrategy.storespecialcodes,
                                sfi->sfxstrategy.spmopt == 0 ? true : false,
                                err);
    if (sfi->bcktab == NULL)
    {
      sfi->leftborder = NULL;
      sfi->numofallcodes = 0;
      haserr = true;
    } else
    {
      sfi->leftborder = gt_bcktab_leftborder(sfi->bcktab);
      sfi->numofallcodes = gt_bcktab_numofallcodes(sfi->bcktab);
      if (sfi->sfxstrategy.spmopt > 0 && prefixlength > 1U)
      {
        GT_INITBITTAB(sfi->markwholeleafbuckets,sfi->numofallcodes);
        estimatedspace += sizeof (*sfi->markwholeleafbuckets *
                                  GT_NUMOFINTSFORBITS(sfi->numofallcodes));
      }
      sizeofbcktab = gt_bcktab_sizeoftable(sfi->numofchars,prefixlength,
                                           sfi->totallength+1,
                                           sfi->sfxstrategy.spmopt == 0
                                             ? true : false);
      estimatedspace += (size_t) sizeofbcktab +
                        gt_bcktab_sizeofworkspace(prefixlength);
      gt_logger_log(logger,"sizeof(bcktab)=" Formatuint64_t,sizeofbcktab);
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
        if (sfi->sfxstrategy.spmopt == 0)
        {
          updateleftborder_getencseqkmers_twobitencoding(encseq,
                                                         readmode,
                                                         prefixlength,
                                                         prefixlength,
                                                         sfi,
                                                         sfi);
        } else
        {
          spmopt_updateleftborder_getencseqkmers_twobitencoding(
                                                        encseq,
                                                        readmode,
                                                        prefixlength,
                                                        sfi->sfxstrategy.spmopt,
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
          sfi->sfxstrategy.spmopt == 0)
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
    gt_bcktab_showleftborder(sfi->bcktab);
#endif
    largestbucketsize
      = gt_bcktab_leftborderpartialsums(&saved_bucketswithoutwholeleaf,
                                        &numofsuffixestosort,
                                        sfi->bcktab,
                                        sfi->markwholeleafbuckets);
    gt_logger_log(sfi->logger, "largest bucket size=%lu",largestbucketsize);

    if (sfi->sfxstrategy.spmopt > 0)
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
      gt_Outlcpinfo_numsuffixes2output_set(sfi->outlcpinfo,
                                           sfi->sfxstrategy.spmopt == 0
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
                                       sfi->totallength,
                                       specialcharacters,
                                       numofsuffixestosort,
                                       sfi->sfxstrategy.suftabcompressedbytes,
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
                                         sfi->bcktab,
                                         numofsuffixestosort,
                                         specialcharacters + 1,
                                         logger);
    gt_assert(sfi->suftabparts != NULL);
    sfi->suffixsortspace
      = gt_suffixsortspace_new(stpgetlargestwidth(sfi->suftabparts),
                               sfi->totallength,
                               sfi->sfxstrategy.suftabcompressedbytes);
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
    if (numofparts > 1U)
    {
      FILE *bcktmpfilefp;

      gt_assert(sfi != NULL);
      sfi->bcktmpfilename = gt_str_new();
      bcktmpfilefp = gt_xtmpfp(sfi->bcktmpfilename);
      gt_logger_log(logger, "store bcktab in \"%s\"",
                    gt_str_get(sfi->bcktmpfilename));
      if (gt_bcktab_flush_to_file(bcktmpfilefp,sfi->bcktab,err) != 0)
      {
        haserr = true;
      }
      gt_fa_fclose(bcktmpfilefp);
    } else
    {
      gt_assert(sfi != NULL);
      sfi->bcktmpfilename = NULL;
    }
  }
  if (haserr)
  {
    (void) gt_Sfxiterator_delete(sfi,NULL);
    return NULL;
  }
  return sfi;
}

static void preparethispart(Sfxiterator *sfi)
{
  unsigned long partwidth;
  unsigned int numofparts = stpgetnumofparts(sfi->suftabparts);
  GtBucketspec2 *bucketspec2 = NULL;

  if (sfi->part == 0 && sfi->withprogressbar)
  {
    gt_progressbar_start(&sfi->bucketiterstep,
                         (unsigned long long) sfi->numofallcodes);
  }
  sfi->currentmincode = stpgetcurrentmincode(sfi->part,sfi->suftabparts);
  sfi->currentmaxcode = stpgetcurrentmaxcode(sfi->part,sfi->suftabparts);
  sfi->widthofpart = stpgetcurrentwidthofpart(sfi->part,sfi->suftabparts);
  if (stpgetnumofparts(sfi->suftabparts) > 1U)
  {
    gt_logger_log(sfi->logger,"compute part %u: "
                              "%lu suffixes,%lu buckets from "
                              "%lu(%lu)..%lu(%lu)",
                  sfi->part,
                  stpgetcurrentwidthofpart(sfi->part,sfi->suftabparts),
                  sfi->currentmaxcode - sfi->currentmincode + 1,
                  sfi->currentmincode,
                  sfi->currentmincode * gt_bcktab_sizeofbasetype(sfi->bcktab),
                  sfi->currentmaxcode,
                  sfi->currentmaxcode * gt_bcktab_sizeofbasetype(sfi->bcktab));
    gt_bcktab_assignboundsforpart(sfi->bcktab,
                                  gt_str_get(sfi->bcktmpfilename),
                                  sfi->part,
                                  sfi->currentmincode,
                                  sfi->currentmaxcode);
  }
  gt_suffixsortspace_offset_set(sfi->suffixsortspace,
                                stpgetcurrentsuftaboffset(sfi->part,
                                                          sfi->suftabparts));
  if (sfi->sfxprogress != NULL)
  {
    gt_timer_show_progress(sfi->sfxprogress, "inserting suffixes into buckets",
                           stdout);
  }
  if (sfi->sfxstrategy.spmopt == 0)
  {
    if (sfi->sfxstrategy.storespecialcodes)
    {
      sfx_derivespecialcodesfromtable(sfi,(numofparts == 1U) ? true : false);
    } else
    {
      sfx_derivespecialcodesonthefly(sfi);
    }
  }
  if (sfi->prefixlength > 1U
      && gt_has_twobitencoding(sfi->encseq)
      && !sfi->sfxstrategy.kmerswithencseqreader)
  {
    insertsuffix_getencseqkmers_twobitencoding(
                                           sfi->encseq,
                                           sfi->readmode,
                                           sfi->prefixlength,
                                           sfi->sfxstrategy.spmopt == 0
                                             ? sfi->prefixlength
                                             : sfi->sfxstrategy.spmopt,
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
  if (sfi->sfxprogress != NULL)
  {
    gt_timer_show_progress(sfi->sfxprogress, "sorting the buckets", stdout);
  }
  /* exit(0); just for testing */
  partwidth = stpgetcurrentsumofwdith(sfi->part,sfi->suftabparts);

  if (numofparts == 1U && sfi->outlcpinfo == NULL &&
      sfi->prefixlength >= 2U && sfi->sfxstrategy.spmopt == 0)
  {
    bucketspec2 = gt_copysort_new(sfi->bcktab,sfi->encseq,sfi->readmode,
                                  partwidth,sfi->numofchars);
  }
  if (sfi->sfxstrategy.differencecover > 0)
  {
    gt_differencecoversetsuffixsortspace(sfi->dcov,sfi->suffixsortspace);
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
    preparethispart(sfi);
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
  if (sfi->sfxstrategy.spmopt == 0)
  {
    gt_suffixsortspace_offset_set(sfi->suffixsortspace,
                                  sfi->totallength - sfi->specialcharacters);
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

int gt_Sfxiterator_bcktab2file(FILE *fp, const Sfxiterator *sfi, GtError *err)
{
  gt_error_check(err);
  return gt_bcktab_flush_to_file(fp,sfi->bcktab,err);
}

unsigned long gt_Sfxiterator_longest(const Sfxiterator *sfi)
{
  gt_assert(sfi != NULL);
  if (sfi->sfxstrategy.spmopt == 0)
  {
    return gt_suffixsortspace_longest(sfi->suffixsortspace);
  } else
  {
    return 0;
  }
}
