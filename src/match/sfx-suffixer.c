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
#include "intcode-def.h"
#include "spacedef.h"
#include "esa-fileend.h"
#include "suffixptr.h"
#include "sfx-diffcov.h"
#include "sfx-partssuf.h"
#include "sfx-suffixer.h"
#include "sfx-enumcodes.h"
#include "sfx-strategy.h"
#include "sfx-copysort.h"
#include "sfx-mappedstr.h"
#include "sfx-bentsedg.h"
#include "stamp.h"

GT_DECLAREARRAYSTRUCT(Suffixptr);

struct Sfxiterator
{
  bool storespecials;
  GtCodetype currentmincode,
           currentmaxcode;
  unsigned long specialcharacters,
         widthofpart,
         totallength;
  Suffixsortspace suffixsortspace;
  Definedunsignedlong longest;
  unsigned long nextfreeCodeatposition;
  Codeatposition *spaceCodeatposition;
  Suftabparts *suftabparts;
  const GtEncseq *encseq;
  GtReadmode readmode;
  Outlcpinfo *outlcpinfo;
  unsigned int part,
               numofchars,
               prefixlength;
  GtArraySuffixptr fusp;
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

    ecp = gt_newEnumcodeatposition(encseq,
                                readmode,
                                prefixlength,
                                numofchars);
    for (insertindex = 0; gt_nextEnumcodeatposition(&specialcontext,ecp);
         insertindex++)
    {
      codelist[insertindex].maxprefixindex = specialcontext.maxprefixindex;
      codelist[insertindex].position = specialcontext.position;
      codelist[insertindex].code
        = gt_computefilledqgramcode(ecp,
                                 specialcontext.maxprefixindex,
                                 specialcontext.position -
                                 specialcontext.maxprefixindex);
    }
    gt_freeEnumcodeatposition(&ecp);
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
    fprintf(stderr,"len1 = %lu != %lu = len2\n",len1,len2);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  for (idx=0; idx<len1; idx++)
  {
    if (codelist1[idx].position != codelist2[idx].position)
    {
      fprintf(stderr,"idx %lu, codelist1.position = " FormatSeqpos " != "
                      FormatSeqpos " = codelist2.position\n",idx,
                      PRINTSeqposcast(codelist1[idx].position),
                      PRINTSeqposcast(codelist2[idx].position));
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

  ALLOCASSIGNSPACE(spaceCodeatposition2,NULL,Codeatposition,
                   realspecialranges+1);
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
  FREESPACE(spaceCodeatposition2);
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
    cc = gt_encseq_get_encoded_char_nospecial(encseq,pos + idx,
                                                      readmode);
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

  if (kmercode->definedspecialpos)
  {
    if (sfi->storespecials)
    {
      if (kmercode->specialpos > 0)
      {
        if (sfi->sfxstrategy.storespecialcodes)
        {
          Codeatposition *cp;

          cp = sfi->spaceCodeatposition + sfi->nextfreeCodeatposition++;
          gt_assert(kmercode->code <= (GtCodetype) MAXCODEVALUE);
          cp->code = (unsigned int) kmercode->code;
          gt_assert(kmercode->specialpos <= MAXPREFIXLENGTH);
          cp->maxprefixindex = kmercode->specialpos;
          cp->position = position + kmercode->specialpos;
        }
        sfi->storespecials = false;
        gt_assert(kmercode->code > 0);
        sfi->leftborder[kmercode->code]++;
      }
    } else
    {
      if (kmercode->specialpos > 0)
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
        printf("### position " FormatSeqpos ", code2 = %lu != 0\n",
                        PRINTSeqposcast(position),code2);
        printf("previouscode = " FormatGtCodetype "\n",
                        previouscode);
        if (previouskmercodedefined)
        {
          printf("previouskmercodedefined = true\n");
          printf("previousstorespecials = %s\n",
                  previousstorespecials ? "true" : "false");
          printf("previousspecialpos = %u\n",previousspecialpos);
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
  previousspecialpos = kmercode->specialpos;
#endif
}

static void insertwithoutspecial(void *processinfo,
                                 unsigned long position,
                                 const GtKmercode *kmercode)
{
  if (!kmercode->definedspecialpos)
  {
    Sfxiterator *sfi = (Sfxiterator *) processinfo;

    if (kmercode->code >= sfi->currentmincode &&
        kmercode->code <= sfi->currentmaxcode)
    {
      unsigned long stidx = --sfi->leftborder[kmercode->code];
      suffixptrset2(&sfi->suffixsortspace,stidx,position);
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
          suffixptrset2(&sfi->suffixsortspace,stidx,
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
    ecp = gt_newEnumcodeatposition(sfi->encseq,sfi->readmode,
                                sfi->prefixlength,
                                sfi->numofchars);
    while (gt_nextEnumcodeatposition(&specialcontext,ecp))
    {
      if (prefixindex <= specialcontext.maxprefixindex)
      {
        if (gt_computefilledqgramcodestopatmax(&code,
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
            suffixptrset2(&sfi->suffixsortspace,stidx,
                          specialcontext.position - prefixindex);
          }
        }
      }
    }
    gt_freeEnumcodeatposition(&ecp);
  }
}

void gt_freeSfxiterator(Sfxiterator **sfiptr)
{
  Sfxiterator *sfi = (Sfxiterator *) *sfiptr;
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
  FREESPACE(sfi->spaceCodeatposition);
  FREESPACE(sfi->suffixsortspace.sortspace);
  gt_freesuftabparts(sfi->suftabparts);
  if (sfi->bcktab != NULL)
  {
    gt_bcktab_delete(&sfi->bcktab);
  }
  if (sfi->dcov != NULL)
  {
    gt_differencecover_delete(sfi->dcov);
  }
  FREESPACE(*sfiptr);
}

#ifdef SKDEBUG
static void showleftborder(const unsigned long *leftborder,
                           GtCodetype numofallcodes)
{
  GtCodetype i;

  for (i=0; i<MIN(numofallcodes,(GtCodetype) 1024); i++)
  {
    printf("leftborder[" FormatGtCodetype "]=" FormatSeqpos "\n",
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

  kmercodeiterator = gt_kmercodeiterator_encseq_new(encseq,readmode,kmersize);
  if (!gt_kmercodeiterator_inputexhausted(kmercodeiterator))
  {
    unsigned long position = 0;

    while ((kmercodeptr = gt_kmercodeiterator_encseq_next(kmercodeiterator))
                          != NULL)
    {
      updatekmercount(sfi,position++,kmercodeptr);
    }
    gt_kmercodeiterator_delete(kmercodeiterator);
  }
}

void getencseqkmersinsertwithoutspecial(const GtEncseq *encseq,
                                        GtReadmode readmode,
                                        unsigned int kmersize,
                                        Sfxiterator *sfi)
{
  GtKmercodeiterator *kmercodeiterator;
  const GtKmercode *kmercodeptr;

  kmercodeiterator = gt_kmercodeiterator_encseq_new(encseq,readmode,kmersize);
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

Sfxiterator *gt_newSfxiterator(const GtEncseq *encseq,
                               GtReadmode readmode,
                               unsigned int prefixlength,
                               unsigned int numofparts,
                               void *voidoutlcpinfo,
                               const Sfxstrategy *sfxstrategy,
                               GtProgressTimer *sfxprogress,
                               bool withprogressbar,
                               GtLogger *logger,
                               GtError *err)
{
  Sfxiterator *sfi = NULL;
  unsigned long realspecialranges, specialcharacters;
  bool haserr = false;

  gt_error_check(err);

  realspecialranges = gt_encseq_realspecialranges(encseq);
  specialcharacters = gt_encseq_specialcharacters(encseq);
  gt_assert(prefixlength > 0);
  if (sfxstrategy != NULL && sfxstrategy->storespecialcodes &&
      prefixlength > MAXPREFIXLENGTH)
  {
    gt_error_set(err,"argument for option -pl must be in the range [1,%u]",
                  MAXPREFIXLENGTH);
    haserr = true;
  }
  if (!haserr)
  {
    ALLOCASSIGNSPACE(sfi,NULL,Sfxiterator,1);
    if (sfxstrategy != NULL && sfxstrategy->storespecialcodes)
    {
      ALLOCASSIGNSPACE(sfi->spaceCodeatposition,NULL,
                       Codeatposition,realspecialranges+1);
      gt_logger_log(logger,"sizeof (spaceCodeatposition)=%lu",
                              (unsigned long) (sizeof (Codeatposition) *
                                               (realspecialranges+1)));
    } else
    {
      sfi->spaceCodeatposition = NULL;
    }
    sfi->bcktab = NULL;
    sfi->nextfreeCodeatposition = 0;
    sfi->suffixsortspace.sortspace = NULL;
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
       if (sfxstrategy->cmpcharbychar
             || !gt_encseq_bitwise_cmp_ok(encseq))
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
    gt_logger_log(logger,"totallength=%lu",
                        sfi->totallength);
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
                              logger,
                              err);
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
    gt_assert(sfi != NULL);
    sfi->storespecials = true;
    if (sfxprogress != NULL)
    {
      gt_progress_timer_start_new_state(sfxprogress,
                                        "counting prefix distribution",
                                        stdout);
    }
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
      reversespecialcodes(sfi->spaceCodeatposition,sfi->nextfreeCodeatposition);
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
    gt_assert(sfi->leftborder != NULL);
#ifdef SKDEBUG
    showleftborder(sfi->leftborder,sfi->numofallcodes);
#endif
    gt_bcktab_leftborderpartialsums(sfi->bcktab,
                                 sfi->totallength - specialcharacters);
    sfi->suftabparts = gt_newsuftabparts(numofparts,
                                         sfi->leftborder,
                                         sfi->numofallcodes,
                                         sfi->totallength - specialcharacters,
                                         specialcharacters + 1,
                                         logger);
    gt_assert(sfi->suftabparts != NULL);
    ALLOCASSIGNSPACE(sfi->suffixsortspace.sortspace,NULL,GtUlong,
                     stpgetlargestwidth(sfi->suftabparts));
    sfi->longest.defined = false;
    sfi->longest.valueunsignedlong = 0;
    if (gt_encseq_has_specialranges(sfi->encseq))
    {
      sfi->sri = gt_specialrangeiterator_new(sfi->encseq,
                                         GT_ISDIRREVERSE(sfi->readmode)
                                           ? false : true);
    } else
    {
      sfi->sri = NULL;
    }
    sfi->fusp.spaceSuffixptr = sfi->suffixsortspace.sortspace;
    sfi->fusp.allocatedSuffixptr = stpgetlargestwidth(sfi->suftabparts);
    sfi->overhang.start = sfi->overhang.end = 0;
  }
  if (haserr)
  {
    if (sfi != NULL)
    {
      gt_freeSfxiterator(&sfi);
    }
    return NULL;
  }
  return sfi;
}

bool gt_sfi2longestsuffixpos(unsigned long *longest,const Sfxiterator *sfi)
{
  if (sfi->longest.defined)
  {
    *longest = sfi->longest.valueunsignedlong;
    return true;
  }
  return false;
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
  sfi->suffixsortspace.sortspaceoffset
    = stpgetcurrentsuftaboffset(sfi->part,sfi->suftabparts);
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
      gt_assert(sfi->suffixsortspace.sortspace != NULL);
      gt_qsufsort(sfi->suffixsortspace.sortspace,
                  partwidth,
                  -1,
                  NULL,
                  &sfi->longest.valueunsignedlong,
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
      sfi->longest.defined = true;
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
                              sfi->suffixsortspace.sortspace -
                              sfi->suffixsortspace.sortspaceoffset,
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
                              dc_sortunsortedbucket,
                              sfi->logger);
    } else
    {
      gt_sortallbuckets (&sfi->suffixsortspace,
                         &sfi->longest,
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
      gt_assert(sfi->suffixsortspace.sortspaceoffset == 0);
      gt_copysort_derivesorting(bucketspec2,&sfi->suffixsortspace,sfi->logger);
      gt_copysort_delete(bucketspec2);
      bucketspec2 = NULL;
    }
  }
  sfi->part++;
}

int gt_postsortsuffixesfromstream(Sfxiterator *sfi, const GtStr *indexname,
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
  FREESPACE(sfi->suffixsortspace.sortspace);
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
                &sfi->longest.valueunsignedlong,
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
    sfi->longest.defined = true;
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
      SUFFIXPTRSET(sfi->fusp.spaceSuffixptr,sfi->fusp.nextfreeSuffixptr,
                   GT_REVERSEPOS(sfi->totallength,pos));
      sfi->fusp.nextfreeSuffixptr++;
      if (pos == leftpos)
      {
        break;
      }
      pos--;
    } else
    {
      SUFFIXPTRSET(sfi->fusp.spaceSuffixptr,sfi->fusp.nextfreeSuffixptr,pos);
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
          SUFFIXPTRSET(sfi->fusp.spaceSuffixptr,sfi->fusp.nextfreeSuffixptr,
                       sfi->totallength);
          sfi->fusp.nextfreeSuffixptr++;
          sfi->exhausted = true;
        }
        break;
      }
    }
  }
}

const void *gt_nextSfxiterator(unsigned long *numberofsuffixes, /* XXX */
                               bool *specialsuffixes,
                               Sfxiterator *sfi)
{
  if (sfi->part < stpgetnumofparts(sfi->suftabparts))
  {
    preparethispart(sfi);
    *numberofsuffixes = sfi->widthofpart;
    *specialsuffixes = false;
    return sfi->suffixsortspace.sortspace;
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
  fillspecialnextpage(sfi);
  gt_assert(sfi->fusp.nextfreeSuffixptr > 0);
  *numberofsuffixes = sfi->fusp.nextfreeSuffixptr;
  *specialsuffixes = true;
  return (void *) sfi->suffixsortspace.sortspace;
}

int gt_sfibcktab2file(FILE *fp,
                   const Sfxiterator *sfi,
                   GtError *err)
{
  gt_error_check(err);
  return gt_bcktab2file(fp,sfi->bcktab,err);
}

unsigned int getprefixlenbits(void)
{
  return (unsigned int) PREFIXLENBITS;
}
