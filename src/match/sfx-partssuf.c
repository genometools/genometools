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

#include "core/assert_api.h"
#include "core/ma_api.h"
#include "core/codetype.h"
#include "core/spacecalc.h"
#include "core/format64.h"
#include "core/log.h"
#include "bcktab.h"
#include "sfx-suffixgetset.h"
#include "spmsuftab.h"
#include "sfx-partssuf.h"
#include "firstcodes-tab.h"

typedef struct
{
  unsigned long nextidx,
                widthofpart,
                suftaboffset,
                sumofwidth,
                minindex,
                maxindex;
} GtSuftabpartcomponent;

struct GtSuftabparts
{
  GtSuftabpartcomponent *components;
  bool indexrange_available;
  unsigned int numofparts;
  unsigned long largestsizemappedpartwise,
                largestsuftabwidth;
  const GtFirstcodestab *fct;
};

void gt_suftabparts_showallrecords(const GtSuftabparts *suftabparts,
                                   bool withminmaxindex)
{
  unsigned int part;
  unsigned long totalwidth;

  gt_assert(suftabparts != NULL);
  gt_assert(suftabparts->numofparts >= 1U);
  totalwidth = suftabparts->components[suftabparts->numofparts - 1].sumofwidth;
  for (part = 0; part < suftabparts->numofparts; part++)
  {
    if (withminmaxindex)
    {
      gt_log_log("part %u: width=%lu (%.2f%%) offset=%lu nextidx=%lu "
                 "minindex=%lu maxindex=%lu ",
                 part,
                 suftabparts->components[part].widthofpart,
                 100.00 * (double) suftabparts->components[part].widthofpart/
                                   totalwidth,
                 suftabparts->components[part].suftaboffset,
                 suftabparts->components[part].nextidx,
                 gt_suftabparts_minindex(part,suftabparts),
                 gt_suftabparts_maxindex(part,suftabparts));
    } else
    {
      gt_log_log("part %u: width=%lu (%.2f%%) offset=%lu nextidx=%lu",
                 part,
                 suftabparts->components[part].widthofpart,
                 100.00 * (double) suftabparts->components[part].widthofpart/
                                   totalwidth,
                 suftabparts->components[part].nextidx,
                 suftabparts->components[part].suftaboffset);
    }
  }
  gt_log_log("variance %.0f",gt_suftabparts_variance(suftabparts));
}

static GtCodetype gt_suftabparts_minindex_raw(unsigned int part,
                                              const GtSuftabparts *suftabparts)
{
  GtCodetype minindex;
  gt_assert(suftabparts != NULL && part < suftabparts->numofparts);

  minindex = (part == 0) ? 0 : suftabparts->components[part-1].nextidx + 1;
  if (suftabparts->fct != NULL)
  {
    return gt_firstcodes_sample2full(suftabparts->fct,minindex);
  }
  return minindex;
}

static GtCodetype gt_suftabparts_maxindex_raw(unsigned int part,
                                              const GtSuftabparts *suftabparts)
{
  GtCodetype maxindex;

  gt_assert(suftabparts != NULL && part < suftabparts->numofparts);
  if (suftabparts->fct != NULL)
  {
    maxindex = suftabparts->components[part].nextidx;
    return gt_firstcodes_sample2full(suftabparts->fct,maxindex);
  }
  maxindex = (part == suftabparts->numofparts - 1)
               ? suftabparts->components[part].nextidx - 1
               : suftabparts->components[part].nextidx;
  return maxindex;
}

static void gt_suftabparts_removeemptyparts(GtSuftabparts *suftabparts,
                                            GT_UNUSED unsigned long totalwidth,
                                            GtLogger *logger)
{
#undef SKDEBUG
#ifdef SKDEBUG
  gt_log_log("before removal");
  gt_suftabparts_showallrecords(suftabparts,false);
#endif
  gt_assert(suftabparts != NULL);
  if (suftabparts->numofparts > 0)
  {
    unsigned int destpart, srcpart;
    unsigned long sumwidth = 0;

    for (destpart = 0, srcpart = 0; srcpart < suftabparts->numofparts;
         srcpart++)
    {
      if (suftabparts->components[srcpart].widthofpart > 0)
      {
        if (destpart < srcpart)
        {
          suftabparts->components[destpart] = suftabparts->components[srcpart];
        }
        destpart++;
      }
    }
    if (destpart < srcpart)
    {
      suftabparts->components[destpart-1].nextidx
        = suftabparts->components[suftabparts->numofparts-1].nextidx;
      suftabparts->numofparts -= (srcpart - destpart);
      gt_assert(suftabparts->numofparts == destpart);
    }
    for (srcpart = 0; srcpart < suftabparts->numofparts; srcpart++)
    {
      gt_assert(suftabparts->components[srcpart].widthofpart > 0);
      sumwidth+=suftabparts->components[srcpart].widthofpart;
      gt_logger_log(logger,"widthofpart[%u]=%lu",
                    srcpart,
                    suftabparts->components[srcpart].widthofpart);
    }
    gt_assert(sumwidth == totalwidth);
  }
#ifdef SKDEBUG
  gt_log_log("after removal");
  gt_suftabparts_showallrecords(suftabparts,true);
#endif
}

GtSuftabparts *gt_suftabparts_new(unsigned int numofparts,
                                  const GtBcktab *bcktab,
                                  const GtFirstcodestab *fct,
                                  const GtSfxmappedrangelist *sfxmrlist,
                                  unsigned long numofsuffixestoinsert,
                                  unsigned long fullspecials,
                                  GtLogger *logger)
{
  GtSuftabparts *suftabparts;
  unsigned long size_mapped;
  unsigned int part;

  gt_assert((bcktab == NULL && fct != NULL) ||
            (bcktab != NULL && fct == NULL));
  suftabparts = gt_malloc(sizeof *suftabparts);
  suftabparts->largestsizemappedpartwise = 0;
  suftabparts->fct = fct;
  gt_assert(suftabparts != NULL);
  if (numofsuffixestoinsert == 0)
  {
    suftabparts->numofparts = 0;
  } else
  {
    if (numofsuffixestoinsert < (unsigned long) numofparts ||
        (bcktab != NULL && gt_bcktab_prefixlength(bcktab) == 1U))
    {
      suftabparts->numofparts = 1U;
    } else
    {
      suftabparts->numofparts = numofparts;
    }
  }
  if (suftabparts->numofparts == 0)
  {
    suftabparts->largestsuftabwidth = fullspecials/numofparts+1;
    suftabparts->largestsizemappedpartwise
      = gt_Sfxmappedrangelist_size_entire(sfxmrlist);
    suftabparts->components = NULL;
  } else
  {
    unsigned int remainder;
    unsigned long secondidx, firstbound = 0, secondbound,
                  suftaboffset = 0, sumofwidth = 0;
    const unsigned long widthofsuftabpart
      = numofsuffixestoinsert/suftabparts->numofparts;

    suftabparts->components
      = gt_malloc(sizeof (*suftabparts->components) * suftabparts->numofparts);
    remainder = (unsigned int) (numofsuffixestoinsert %
                                (unsigned long) suftabparts->numofparts);
    suftabparts->largestsuftabwidth = 0;
    for (part=0; part < suftabparts->numofparts; part++)
    {
      if (remainder > 0)
      {
        suftaboffset += widthofsuftabpart + 1;
        remainder--;
      } else
      {
        suftaboffset += widthofsuftabpart;
      }
      if (part == suftabparts->numofparts - 1)
      {
        secondidx = bcktab != NULL ? gt_bcktab_numofallcodes(bcktab)
                                   : gt_firstcodes_numofsamples(fct);
      } else
      {
        secondidx = bcktab != NULL
                      ? gt_bcktab_findfirstlarger(bcktab,suftaboffset)
                      : gt_firstcodes_findfirstsamplelarger(fct,suftaboffset);
      }
      suftabparts->components[part].nextidx = secondidx;
      secondbound = bcktab != NULL
                      ? gt_bcktab_get_leftborder(bcktab,secondidx)
                      : gt_firstcodes_get_sample(fct,secondidx);
      if (part == 0)
      {
        suftabparts->components[part].widthofpart = secondbound;
        suftabparts->components[part].suftaboffset = 0;
      } else
      {
        gt_assert(secondbound >= firstbound);
        suftabparts->components[part].widthofpart = secondbound - firstbound;
        suftabparts->components[part].suftaboffset = firstbound;
      }
      if (suftabparts->largestsuftabwidth <
          suftabparts->components[part].widthofpart)
      {
        suftabparts->largestsuftabwidth
          = suftabparts->components[part].widthofpart;
      }
      sumofwidth += suftabparts->components[part].widthofpart;
      suftabparts->components[part].sumofwidth = sumofwidth;
      firstbound = secondbound;
    }
    gt_assert(sumofwidth == numofsuffixestoinsert);
  }
  gt_suftabparts_removeemptyparts(suftabparts,numofsuffixestoinsert,logger);
  if (suftabparts->numofparts > 0)
  {
    gt_assert(suftabparts->components != NULL);
    for (part=0; part < suftabparts->numofparts; part++)
    {
      suftabparts->components[part].minindex
        = gt_suftabparts_minindex_raw(part,suftabparts);
      suftabparts->components[part].maxindex
        = gt_suftabparts_maxindex_raw(part,suftabparts);
    }
    for (part=1U; part < suftabparts->numofparts; part++)
    {
      if (suftabparts->components[part].minindex !=
            suftabparts->components[part-1].maxindex + 1)
      {
        suftabparts->components[part].minindex
          = suftabparts->components[part-1].maxindex + 1;
      }
    }
    for (part=0; part < suftabparts->numofparts; part++)
    {
      size_mapped = gt_Sfxmappedrangelist_size_mapped(sfxmrlist,
                                    gt_suftabparts_minindex(part,suftabparts),
                                    gt_suftabparts_maxindex(part,suftabparts));
      if (suftabparts->largestsizemappedpartwise < size_mapped)
      {
        suftabparts->largestsizemappedpartwise = size_mapped;
      }
    }
  }
  return suftabparts;
}

double gt_suftabparts_variance(const GtSuftabparts *suftabparts)
{

  gt_assert(suftabparts->numofparts > 0);
  if (suftabparts->numofparts == 1U)
  {
    return 0.0;
  } else
  {
    double meanwidth, difference, sum = 0.0;
    unsigned int part;

    meanwidth = (double)
                suftabparts->components[suftabparts->numofparts-1].sumofwidth/
                suftabparts->numofparts;
    for (part = 0; part < suftabparts->numofparts; part++)
    {
      difference = (double) suftabparts->components[part].widthofpart
                   - meanwidth;
      sum += difference * difference;
    }
    return sum/(double) (suftabparts->numofparts-1);
  }
}

void gt_suftabparts_delete(GtSuftabparts *suftabparts)
{
  if (suftabparts != NULL)
  {
    gt_free(suftabparts->components);
    gt_free(suftabparts);
  }
}

GtCodetype gt_suftabparts_minindex(unsigned int part,
                                   const GtSuftabparts *suftabparts)
{
  return suftabparts->components[part].minindex;
}

GtCodetype gt_suftabparts_maxindex(unsigned int part,
                                   const GtSuftabparts *suftabparts)
{
  return suftabparts->components[part].maxindex;
}

GtCodetype gt_suftabparts_maxindex_last(const GtSuftabparts *suftabparts)
{
  gt_assert(suftabparts->numofparts > 0);
  return gt_suftabparts_maxindex(suftabparts->numofparts-1,suftabparts);
}

unsigned long gt_suftabparts_offset(unsigned int part,
                                    const GtSuftabparts *suftabparts)
{
  gt_assert(suftabparts != NULL && part < suftabparts->numofparts);
  return suftabparts->components[part].suftaboffset;
}

unsigned long gt_suftabparts_sumofwidth(unsigned int part,
                                        const GtSuftabparts *suftabparts)
{
  gt_assert(suftabparts != NULL && part < suftabparts->numofparts);
  return suftabparts->components[part].sumofwidth;
}

unsigned long gt_suftabparts_widthofpart(unsigned int part,
                                         const GtSuftabparts *suftabparts)
{
  gt_assert(suftabparts != NULL && part < suftabparts->numofparts);
  return suftabparts->components[part].widthofpart;
}

unsigned long gt_suftabparts_largest_width(const GtSuftabparts *suftabparts)
{
  gt_assert(suftabparts != NULL);
  return suftabparts->largestsuftabwidth;
}

unsigned long gt_suftabparts_largestsizemappedpartwise(
                                          const GtSuftabparts *suftabparts)
{
  gt_assert(suftabparts != NULL);
  return suftabparts->largestsizemappedpartwise;
}

unsigned int gt_suftabparts_numofparts(const GtSuftabparts *suftabparts)
{
  gt_assert(suftabparts != NULL);
  return suftabparts->numofparts;
}

int gt_suftabparts_fit_memlimit(size_t estimatedspace,
                                unsigned long maximumspace,
                                const GtBcktab *bcktab,
                                const GtFirstcodestab *fct,
                                const GtSfxmappedrangelist *sfxmrlist,
                                unsigned long totallength,
                                unsigned int bitsforseqnumrelpos,
                                unsigned long specialcharacters,
                                unsigned long numofsuffixestosort,
                                bool suftabuint,
                                GtError *err)
{
  unsigned int parts;
  GtSuftabparts *suftabparts;
  unsigned long size_mapped = gt_Sfxmappedrangelist_size_entire(sfxmrlist);

  gt_error_check(err);
  for (parts = 1U; parts <= 500U; parts++)
  {
    uint64_t suftabsize;
    unsigned long numofentries;

    suftabparts = gt_suftabparts_new(parts,
                                     bcktab,
                                     fct,
                                     sfxmrlist,
                                     numofsuffixestosort,
                                     specialcharacters + 1,
                                     NULL);
    gt_assert(suftabparts != NULL);
    numofentries = gt_suftabparts_largest_width(suftabparts);
    if (bcktab != NULL)
    {
      gt_assert(fct == NULL);
      suftabsize = gt_suffixsortspace_requiredspace(numofentries,
                                                    totallength,
                                                    suftabuint);
    } else
    {
      gt_assert(fct != NULL);
      suftabsize = (uint64_t) gt_spmsuftab_requiredspace(numofentries,
                                                         totallength,
                                                         bitsforseqnumrelpos);
    }
    if (parts == 1U)
    {
      if (suftabsize + (uint64_t) estimatedspace <= (uint64_t) maximumspace)
      {
        gt_suftabparts_delete(suftabparts);
        return (int) parts;
      }
    } else
    {
      unsigned long largest
        = gt_suftabparts_largestsizemappedpartwise(suftabparts);
      if (suftabsize
          + (uint64_t) largest
          + (uint64_t) estimatedspace
          - (uint64_t) size_mapped <= (uint64_t) maximumspace)
      {
        gt_log_log("return parts = %u as suftabsize=%.2f +"
                   "largest=%.2f + estimated=%.2f - size_mapped=%2.f <= %.2f",
                    parts,
                    GT_MEGABYTES(suftabsize),
                    GT_MEGABYTES(largest),
                    GT_MEGABYTES(estimatedspace),
                    GT_MEGABYTES(size_mapped),
                    GT_MEGABYTES(maximumspace));
        gt_suftabparts_delete(suftabparts);
        return (int) parts;
      }
    }
    gt_suftabparts_delete(suftabparts);
  }
  gt_error_set(err,"cannot compute enhanced suffix array in at most %lu bytes",
                   maximumspace);
  return -1;
}
