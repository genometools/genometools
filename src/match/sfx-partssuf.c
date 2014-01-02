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
  GtUword widthofpart,
          suftaboffset,
          sumofwidth,
          nextidx,
          minindex,
          maxindex;
} GtSuftabpartcomponent;

struct GtSuftabparts
{
  GtSuftabpartcomponent *components;
  bool indexrange_available;
  unsigned int numofparts;
  GtUword codeoffset,
          partoffset,
          largestsizemappedpartwise,
          largestsuftabwidth;
  const GtFirstcodestab *fct;
};

void gt_suftabparts_showallrecords(const GtSuftabparts *suftabparts,
                                   bool withminmaxindex)
{
  gt_assert(suftabparts != NULL);
  if (suftabparts->numofparts >= 1U)
  {
    unsigned int part;
    GtUword totalwidth;

    gt_assert(suftabparts->components != NULL);
    totalwidth = gt_suftabparts_sumofwidth(suftabparts->numofparts - 1,
                                           suftabparts);
    for (part = 0; part < suftabparts->numofparts; part++)
    {
      if (withminmaxindex)
      {
        gt_log_log("part %u: width="GT_WU" (%.2f%%) suftaboffset="GT_WU
                   ", sumofwidth="GT_WU ", minindex="GT_WU" maxindex="GT_WU" ",
                   part,
                   suftabparts->components[part].widthofpart,
                   100.00 * (double) suftabparts->components[part].widthofpart/
                                     totalwidth,
                   suftabparts->components[part].suftaboffset,
                   gt_suftabparts_sumofwidth(part,suftabparts),
                   gt_suftabparts_minindex(part,suftabparts),
                   gt_suftabparts_maxindex(part,suftabparts));
      } else
      {
        gt_log_log("part %u: width="GT_WU" (%.2f%%) suftaboffset="GT_WU" ",
                   part,
                   suftabparts->components[part].widthofpart,
                   100.00 * (double) suftabparts->components[part].widthofpart/
                                     totalwidth,
                   suftabparts->components[part].suftaboffset);
      }
    }
    gt_log_log("variance %.0f",gt_suftabparts_variance(suftabparts));
  }
}

static GtCodetype gt_suftabparts_minindex_raw(unsigned int part,
                                              const GtSuftabparts *suftabparts)
{
  GtCodetype nextidx;

  gt_assert(suftabparts != NULL && part < suftabparts->numofparts);
  nextidx = (part == 0) ? suftabparts->codeoffset
                        : suftabparts->components[part-1].nextidx + 1;
  if (suftabparts->fct != NULL)
  {
    return gt_firstcodes_sample2full(suftabparts->fct,nextidx);
  }
  return nextidx;
}

static GtCodetype gt_suftabparts_maxindex_raw(unsigned int part,
                                              const GtSuftabparts *suftabparts)
{
  GtCodetype nextidx;

  gt_assert(suftabparts != NULL && suftabparts->components != NULL &&
            part < suftabparts->numofparts);
  nextidx = suftabparts->components[part].nextidx;
  if (suftabparts->fct != NULL)
  {
    return gt_firstcodes_sample2full(suftabparts->fct,nextidx);
  }
  return nextidx;
}

static void gt_suftabparts_removeemptyparts(GtSuftabparts *suftabparts,
                                            GT_UNUSED GtUword totalwidth,
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
    GtUword sumwidth = 0;

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
      gt_logger_log(logger,"widthofpart[%u]="GT_WU"",
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
                                  GtCodetype mincode,
                                  GtCodetype maxcode,
                                  const GtFirstcodestab *fct,
                                  const GtSfxmappedrangelist *sfxmrlist,
                                  GtUword numofsuffixestoinsert,
                                  GtUword fullspecials,
                                  GtLogger *logger)
{
  GtSuftabparts *suftabparts;
  GtUword size_mapped;
  unsigned int part;

  gt_assert((bcktab == NULL && fct != NULL) ||
            (bcktab != NULL && fct == NULL));
  suftabparts = gt_malloc(sizeof *suftabparts);
  suftabparts->largestsizemappedpartwise = 0;
  if (mincode <= maxcode && fct == NULL)
  {
    suftabparts->codeoffset = mincode;
    if (mincode > 0)
    {
      suftabparts->partoffset = gt_bcktab_get_leftborder(bcktab,mincode-1);
    } else
    {
      suftabparts->partoffset = 0;
    }
  } else
  {
    suftabparts->codeoffset = 0;
    suftabparts->partoffset = 0;
  }
  suftabparts->fct = fct;
  if (numofsuffixestoinsert == 0)
  {
    suftabparts->numofparts = 0;
  } else
  {
    if (numofsuffixestoinsert < (GtUword) numofparts ||
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
    gt_assert(numofparts != 0);
    suftabparts->largestsuftabwidth = fullspecials/numofparts+1;
    suftabparts->largestsizemappedpartwise
      = gt_Sfxmappedrangelist_size_entire(sfxmrlist);
    suftabparts->components = NULL;
  } else
  {
    unsigned int remainder;
    GtUword firstbound = 0, secondbound, suftaboffset = 0, sumofwidth = 0;
    const GtUword widthofsuftabpart
      = numofsuffixestoinsert/suftabparts->numofparts;

    /*printf("%s: parts = "GT_WUT", mincode = "GT_WU", maxcode = "GT_WU"\n",
             __func__,
             suftabparts->numofparts,
             mincode,maxcode);*/

    suftabparts->components
      = gt_malloc(sizeof (*suftabparts->components) * suftabparts->numofparts);
    remainder = (unsigned int) (numofsuffixestoinsert %
                                (GtUword) suftabparts->numofparts);
    suftabparts->largestsuftabwidth = 0;
    for (part=0; part < suftabparts->numofparts; part++)
    {
      GtUword secondidx;
      if (remainder > 0)
      {
        suftaboffset += widthofsuftabpart + 1;
        remainder--;
      } else
      {
        suftaboffset += widthofsuftabpart;
      }
      if (bcktab != NULL)
      {
        if (part == suftabparts->numofparts - 1)
        {
          if (mincode <= maxcode)
          {
            secondidx = maxcode;
          } else
          {
            secondidx = gt_bcktab_numofallcodes(bcktab) - 1;
          }
        } else
        {
          secondidx = gt_bcktab_findfirstlarger(bcktab,mincode,maxcode,
                                                suftaboffset);
        }
        secondbound = gt_bcktab_get_leftborder(bcktab,secondidx);
        gt_assert(secondbound >= suftabparts->partoffset);
        secondbound -= suftabparts->partoffset;
      } else
      {
        gt_assert(fct != NULL && mincode > maxcode);
        if (part == suftabparts->numofparts - 1)
        {
          secondidx = gt_firstcodes_numofsamples(fct);
        } else
        {
          secondidx = gt_firstcodes_findfirstsamplelarger(fct,suftaboffset);
        }
        secondbound = gt_firstcodes_get_sample(fct,secondidx);
      }
      suftabparts->components[part].nextidx = secondidx;
      if (part == 0)
      {
        suftabparts->components[part].widthofpart = secondbound;
        suftabparts->components[part].suftaboffset = suftabparts->partoffset;
      } else
      {
        gt_assert(secondbound >= firstbound);
        suftabparts->components[part].widthofpart = secondbound - firstbound;
        suftabparts->components[part].suftaboffset
          = firstbound + suftabparts->partoffset;
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
    if (sumofwidth != numofsuffixestoinsert)
    {
      printf("sumofwidth = " GT_WU " != " GT_WU "= numofsuffixestosort\n",
              sumofwidth,numofsuffixestoinsert);
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

    meanwidth
      = (double)
        gt_suftabparts_sumofwidth(suftabparts->numofparts-1,suftabparts)/
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
  gt_assert(suftabparts != NULL && suftabparts->components != NULL &&
            part < suftabparts->numofparts);
  return suftabparts->components[part].minindex;
}

GtCodetype gt_suftabparts_maxindex(unsigned int part,
                                   const GtSuftabparts *suftabparts)
{
  gt_assert(suftabparts != NULL && suftabparts->components != NULL &&
            part < suftabparts->numofparts);
  return suftabparts->components[part].maxindex;
}

GtCodetype gt_suftabparts_maxindex_last(const GtSuftabparts *suftabparts)
{
  gt_assert(suftabparts != NULL && suftabparts->numofparts > 0);
  return gt_suftabparts_maxindex(suftabparts->numofparts-1,suftabparts);
}

GtUword gt_suftabparts_offset(unsigned int part,
                              const GtSuftabparts *suftabparts)
{
  gt_assert(suftabparts != NULL && suftabparts->components != NULL &&
            part < suftabparts->numofparts);
  return suftabparts->components[part].suftaboffset;
}

GtUword gt_suftabparts_sumofwidth(unsigned int part,
                                  const GtSuftabparts *suftabparts)
{
  gt_assert(suftabparts != NULL && suftabparts->components != NULL &&
            part < suftabparts->numofparts);
  return suftabparts->components[part].sumofwidth + suftabparts->partoffset;
}

GtUword gt_suftabparts_widthofpart(unsigned int part,
                                   const GtSuftabparts *suftabparts)
{
  gt_assert(suftabparts != NULL && suftabparts->components != NULL &&
            part < suftabparts->numofparts);
  return suftabparts->components[part].widthofpart;
}

GtUword gt_suftabparts_largest_width(const GtSuftabparts *suftabparts)
{
  gt_assert(suftabparts != NULL);
  return suftabparts->largestsuftabwidth;
}

GtUword gt_suftabparts_largestsizemappedpartwise(
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
                                GtUword maximumspace,
                                const GtBcktab *bcktab,
                                const GtFirstcodestab *fct,
                                const GtSfxmappedrangelist *sfxmrlist,
                                GtUword totallength,
                                unsigned int bitsforseqnumrelpos,
                                GtUword specialcharacters,
                                GtUword numofsuffixestosort,
                                bool suftabuint,
                                GtError *err)
{
  unsigned int parts;

  gt_error_check(err);
  for (parts = 1U; parts <= 500U; parts++)
  {
    uint64_t suftabsize;
    GtUword numofentries;
    GtSuftabparts *suftabparts = gt_suftabparts_new(parts,
                                                    bcktab,
                                                    (GtCodetype) 1,
                                                    (GtCodetype) 0,
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
      GtUword largest
        = gt_suftabparts_largestsizemappedpartwise(suftabparts);
      GtUword size_mapped = gt_Sfxmappedrangelist_size_entire(sfxmrlist);

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
  gt_error_set(err,"cannot compute enhanced suffix array in at most "
               GT_WU" bytes", maximumspace);
  return -1;
}
