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
#include "core/log.h"
#include "bcktab.h"
#include "sfx-suffixgetset.h"
#include "spmsuftab.h"
#include "randomcodes-sfx-partssuf.h"
#include "randomcodes-tab.h"

typedef struct
{
  GtUword nextidx,
                widthofpart,
                suftaboffset,
                sumofwidth,
                minindex,
                maxindex;
} GtSuftabpartcomponent_rc;

struct GtSuftabparts_rc
{
  GtSuftabpartcomponent_rc *components;
  bool indexrange_available;
  unsigned int numofparts;
  GtUword largestsizemappedpartwise,
                largestsuftabwidth;
  const GtRandomcodestab *fct;
};

void gt_suftabparts_rc_showallrecords(const GtSuftabparts_rc *suftabparts_rc,
                                   bool withminmaxindex)
{
  unsigned int part;
  GtUword totalwidth;

  gt_assert(suftabparts_rc != NULL);
  gt_assert(suftabparts_rc->numofparts >= 1U);
  totalwidth =
    suftabparts_rc->components[suftabparts_rc->numofparts - 1].sumofwidth;
  for (part = 0; part < suftabparts_rc->numofparts; part++)
  {
    if (withminmaxindex)
    {
      gt_log_log("part %u: width=" GT_WU " (%.2f%%) offset=" GT_WU " nextidx="
                 GT_WU " minindex=" GT_WU " maxindex=" GT_WU " ",
                 part,
                 suftabparts_rc->components[part].widthofpart,
                 100.00 * (double) suftabparts_rc->components[part].widthofpart/
                   totalwidth,
                 suftabparts_rc->components[part].suftaboffset,
                 suftabparts_rc->components[part].nextidx,
                 gt_suftabparts_rc_minindex(part,suftabparts_rc),
                 gt_suftabparts_rc_maxindex(part,suftabparts_rc));
    } else
    {
      gt_log_log("part %u: width=" GT_WU " (%.2f%%) offset=" GT_WU " nextidx="
                 GT_WU,
                 part,
                 suftabparts_rc->components[part].widthofpart,
                 100.00 * (double) suftabparts_rc->components[part].widthofpart/
                   totalwidth,
                 suftabparts_rc->components[part].nextidx,
                 suftabparts_rc->components[part].suftaboffset);
    }
  }
  gt_log_log("variance %.0f",gt_suftabparts_rc_variance(suftabparts_rc));
}

static GtCodetype gt_suftabparts_rc_minindex_raw(
    unsigned int part, const GtSuftabparts_rc *suftabparts_rc)
{
  GtCodetype minindex;
  gt_assert(suftabparts_rc != NULL && part < suftabparts_rc->numofparts);

  minindex = (part == 0) ? 0 : suftabparts_rc->components[part-1].nextidx + 1;
  if (suftabparts_rc->fct != NULL)
  {
    return gt_randomcodes_sample2full(suftabparts_rc->fct,minindex);
  }
  return minindex;
}

static GtCodetype gt_suftabparts_rc_maxindex_raw(
    unsigned int part, const GtSuftabparts_rc *suftabparts_rc)
{
  GtCodetype maxindex;

  gt_assert(suftabparts_rc != NULL && part < suftabparts_rc->numofparts);
  if (suftabparts_rc->fct != NULL)
  {
    maxindex = suftabparts_rc->components[part].nextidx;
    return gt_randomcodes_sample2full(suftabparts_rc->fct,maxindex);
  }
  maxindex = (part == suftabparts_rc->numofparts - 1)
               ? suftabparts_rc->components[part].nextidx - 1
               : suftabparts_rc->components[part].nextidx;
  return maxindex;
}

static void gt_suftabparts_rc_removeemptyparts(GtSuftabparts_rc *suftabparts_rc,
                                            GT_UNUSED GtUword totalwidth,
                                            GtLogger *logger)
{
#undef SKDEBUG
#ifdef SKDEBUG
  gt_log_log("before removal");
  gt_suftabparts_rc_showallrecords(suftabparts_rc,false);
#endif
  gt_assert(suftabparts_rc != NULL);
  if (suftabparts_rc->numofparts > 0)
  {
    unsigned int destpart, srcpart;
    GtUword sumwidth = 0;

    for (destpart = 0, srcpart = 0; srcpart < suftabparts_rc->numofparts;
         srcpart++)
    {
      if (suftabparts_rc->components[srcpart].widthofpart > 0)
      {
        if (destpart < srcpart)
        {
          suftabparts_rc->components[destpart] =
            suftabparts_rc->components[srcpart];
        }
        destpart++;
      }
    }
    if (destpart < srcpart)
    {
      suftabparts_rc->components[destpart-1].nextidx
        = suftabparts_rc->components[suftabparts_rc->numofparts-1].nextidx;
      suftabparts_rc->numofparts -= (srcpart - destpart);
      gt_assert(suftabparts_rc->numofparts == destpart);
    }
    for (srcpart = 0; srcpart < suftabparts_rc->numofparts; srcpart++)
    {
      gt_assert(suftabparts_rc->components[srcpart].widthofpart > 0);
      sumwidth+=suftabparts_rc->components[srcpart].widthofpart;
      gt_logger_log(logger,"widthofpart[%u]=" GT_WU "",
                    srcpart,
                    suftabparts_rc->components[srcpart].widthofpart);
    }
    gt_assert(sumwidth == totalwidth);
  }
#ifdef SKDEBUG
  gt_log_log("after removal");
  gt_suftabparts_rc_showallrecords(suftabparts_rc,true);
#endif
}

GtSuftabparts_rc *gt_suftabparts_rc_new(unsigned int numofparts,
                                  const GtBcktab *bcktab,
                                  const GtRandomcodestab *fct,
                                  const GtSfxmappedrangelist *sfxmrlist,
                                  GtUword numofsuffixestoinsert,
                                  GtUword fullspecials,
                                  GtLogger *logger)
{
  GtSuftabparts_rc *suftabparts_rc;
  GtUword size_mapped;
  unsigned int part;

  gt_assert((bcktab == NULL && fct != NULL) ||
            (bcktab != NULL && fct == NULL));
  suftabparts_rc = gt_malloc(sizeof *suftabparts_rc);
  suftabparts_rc->largestsizemappedpartwise = 0;
  suftabparts_rc->fct = fct;
  gt_assert(suftabparts_rc != NULL);
  if (numofsuffixestoinsert == 0)
  {
    suftabparts_rc->numofparts = 0;
  } else
  {
    if (numofsuffixestoinsert < (GtUword) numofparts ||
        (bcktab != NULL && gt_bcktab_prefixlength(bcktab) == 1U))
    {
      suftabparts_rc->numofparts = 1U;
    } else
    {
      suftabparts_rc->numofparts = numofparts;
    }
  }
  if (suftabparts_rc->numofparts == 0)
  {
    gt_assert(numofparts != 0);
    suftabparts_rc->largestsuftabwidth = fullspecials/numofparts+1;
    suftabparts_rc->largestsizemappedpartwise
      = gt_Sfxmappedrangelist_size_entire(sfxmrlist);
    suftabparts_rc->components = NULL;
  } else
  {
    unsigned int remainder;
    GtUword secondidx, randombound = 0, secondbound,
                  suftaboffset = 0, sumofwidth = 0;
    const GtUword widthofsuftabpart
      = numofsuffixestoinsert/suftabparts_rc->numofparts;

    suftabparts_rc->components
      = gt_malloc(sizeof (*suftabparts_rc->components) *
          suftabparts_rc->numofparts);
    remainder = (unsigned int) (numofsuffixestoinsert %
                                (GtUword) suftabparts_rc->numofparts);
    suftabparts_rc->largestsuftabwidth = 0;
    for (part=0; part < suftabparts_rc->numofparts; part++)
    {
      if (remainder > 0)
      {
        suftaboffset += widthofsuftabpart + 1;
        remainder--;
      } else
      {
        suftaboffset += widthofsuftabpart;
      }
      if (part == suftabparts_rc->numofparts - 1)
      {
        secondidx = bcktab != NULL ? gt_bcktab_numofallcodes(bcktab)
                                   : gt_randomcodes_numofsamples(fct);
      } else
      {
        secondidx = bcktab != NULL
                      ? gt_bcktab_findfirstlarger(bcktab,(GtCodetype) 1,
                                                  (GtCodetype) 0,suftaboffset)
                      : gt_randomcodes_findfirstsamplelarger(fct,suftaboffset);
      }
      suftabparts_rc->components[part].nextidx = secondidx;
      secondbound = bcktab != NULL
                      ? gt_bcktab_get_leftborder(bcktab,secondidx)
                      : gt_randomcodes_get_sample(fct,secondidx);
      if (part == 0)
      {
        suftabparts_rc->components[part].widthofpart = secondbound;
        suftabparts_rc->components[part].suftaboffset = 0;
      } else
      {
        gt_assert(secondbound >= randombound);
        suftabparts_rc->components[part].widthofpart = secondbound -
          randombound;
        suftabparts_rc->components[part].suftaboffset = randombound;
      }
      if (suftabparts_rc->largestsuftabwidth <
          suftabparts_rc->components[part].widthofpart)
      {
        suftabparts_rc->largestsuftabwidth
          = suftabparts_rc->components[part].widthofpart;
      }
      sumofwidth += suftabparts_rc->components[part].widthofpart;
      suftabparts_rc->components[part].sumofwidth = sumofwidth;
      randombound = secondbound;
    }
    gt_assert(sumofwidth == numofsuffixestoinsert);
  }
  gt_suftabparts_rc_removeemptyparts(suftabparts_rc,numofsuffixestoinsert,
      logger);
  if (suftabparts_rc->numofparts > 0)
  {
    gt_assert(suftabparts_rc->components != NULL);
    for (part=0; part < suftabparts_rc->numofparts; part++)
    {
      suftabparts_rc->components[part].minindex
        = gt_suftabparts_rc_minindex_raw(part,suftabparts_rc);
      suftabparts_rc->components[part].maxindex
        = gt_suftabparts_rc_maxindex_raw(part,suftabparts_rc);
    }
    for (part=1U; part < suftabparts_rc->numofparts; part++)
    {
      if (suftabparts_rc->components[part].minindex !=
            suftabparts_rc->components[part-1].maxindex + 1)
      {
        suftabparts_rc->components[part].minindex
          = suftabparts_rc->components[part-1].maxindex + 1;
      }
    }
    for (part=0; part < suftabparts_rc->numofparts; part++)
    {
      size_mapped = gt_Sfxmappedrangelist_size_mapped(sfxmrlist,
          gt_suftabparts_rc_minindex(part,suftabparts_rc),
          gt_suftabparts_rc_maxindex(part,suftabparts_rc));
      if (suftabparts_rc->largestsizemappedpartwise < size_mapped)
      {
        suftabparts_rc->largestsizemappedpartwise = size_mapped;
      }
    }
  }
  return suftabparts_rc;
}

double gt_suftabparts_rc_variance(const GtSuftabparts_rc *suftabparts_rc)
{

  gt_assert(suftabparts_rc->numofparts > 0);
  if (suftabparts_rc->numofparts == 1U)
  {
    return 0.0;
  } else
  {
    double meanwidth, difference, sum = 0.0;
    unsigned int part;

    meanwidth = (double)
      suftabparts_rc->components[suftabparts_rc->numofparts-1].sumofwidth/
      suftabparts_rc->numofparts;
    for (part = 0; part < suftabparts_rc->numofparts; part++)
    {
      difference = (double) suftabparts_rc->components[part].widthofpart
                   - meanwidth;
      sum += difference * difference;
    }
    return sum/(double) (suftabparts_rc->numofparts-1);
  }
}

void gt_suftabparts_rc_delete(GtSuftabparts_rc *suftabparts_rc)
{
  if (suftabparts_rc != NULL)
  {
    gt_free(suftabparts_rc->components);
    gt_free(suftabparts_rc);
  }
}

GtCodetype gt_suftabparts_rc_minindex(unsigned int part,
                                   const GtSuftabparts_rc *suftabparts_rc)
{
  return suftabparts_rc->components[part].minindex;
}

GtCodetype gt_suftabparts_rc_maxindex(unsigned int part,
                                   const GtSuftabparts_rc *suftabparts_rc)
{
  return suftabparts_rc->components[part].maxindex;
}

GtCodetype gt_suftabparts_rc_maxindex_last(
    const GtSuftabparts_rc *suftabparts_rc)
{
  gt_assert(suftabparts_rc->numofparts > 0);
  return gt_suftabparts_rc_maxindex(suftabparts_rc->numofparts-1,
      suftabparts_rc);
}

GtUword gt_suftabparts_rc_offset(unsigned int part,
                                    const GtSuftabparts_rc *suftabparts_rc)
{
  gt_assert(suftabparts_rc != NULL && part < suftabparts_rc->numofparts);
  return suftabparts_rc->components[part].suftaboffset;
}

GtUword gt_suftabparts_rc_sumofwidth(unsigned int part,
                                        const GtSuftabparts_rc *suftabparts_rc)
{
  gt_assert(suftabparts_rc != NULL && part < suftabparts_rc->numofparts);
  return suftabparts_rc->components[part].sumofwidth;
}

GtUword gt_suftabparts_rc_widthofpart(unsigned int part,
                                         const GtSuftabparts_rc *suftabparts_rc)
{
  gt_assert(suftabparts_rc != NULL && part < suftabparts_rc->numofparts);
  return suftabparts_rc->components[part].widthofpart;
}

GtUword gt_suftabparts_rc_largest_width(
    const GtSuftabparts_rc *suftabparts_rc)
{
  gt_assert(suftabparts_rc != NULL);
  return suftabparts_rc->largestsuftabwidth;
}

GtUword gt_suftabparts_rc_largestsizemappedpartwise(
    const GtSuftabparts_rc *suftabparts_rc)
{
  gt_assert(suftabparts_rc != NULL);
  return suftabparts_rc->largestsizemappedpartwise;
}

unsigned int gt_suftabparts_rc_numofparts(
    const GtSuftabparts_rc *suftabparts_rc)
{
  gt_assert(suftabparts_rc != NULL);
  return suftabparts_rc->numofparts;
}

int gt_suftabparts_rc_fit_memlimit(size_t estimatedspace,
                                GtUword maximumspace,
                                const GtBcktab *bcktab,
                                const GtRandomcodestab *fct,
                                const GtSfxmappedrangelist *sfxmrlist,
                                GtUword totallength,
                                unsigned int bitsforseqnumrelpos,
                                GtUword specialcharacters,
                                GtUword numofsuffixestosort,
                                bool suftabuint,
                                GtError *err)
{
  unsigned int parts;
  GtSuftabparts_rc *suftabparts_rc;
  GtUword size_mapped = gt_Sfxmappedrangelist_size_entire(sfxmrlist);

  gt_error_check(err);
  for (parts = 1U; parts <= 500U; parts++)
  {
    size_t suftabsize;
    GtUword numofentries;

    suftabparts_rc = gt_suftabparts_rc_new(parts,
                                     bcktab,
                                     fct,
                                     sfxmrlist,
                                     numofsuffixestosort,
                                     specialcharacters + 1,
                                     NULL);
    gt_assert(suftabparts_rc != NULL);
    numofentries = gt_suftabparts_rc_largest_width(suftabparts_rc);
    if (bcktab != NULL)
    {
      gt_assert(fct == NULL);
      suftabsize = gt_suffixsortspace_requiredspace(numofentries,
                                                    totallength,
                                                    suftabuint);
    } else
    {
      gt_assert(fct != NULL);
      suftabsize = gt_spmsuftab_requiredspace(numofentries,
                                              totallength,
                                              bitsforseqnumrelpos);
    }
    if (parts == 1U)
    {
      if ((GtUword) (suftabsize + estimatedspace) <= maximumspace)
      {
        gt_suftabparts_rc_delete(suftabparts_rc);
        return (int) parts;
      }
    } else
    {
      GtUword largest
        = gt_suftabparts_rc_largestsizemappedpartwise(suftabparts_rc);
      if ((GtUword) (suftabsize + largest + estimatedspace - size_mapped)
                           <= maximumspace)
      {
        gt_log_log("return parts = %u as suftabsize=%.2f +"
                   "largest=%.2f + estimated=%.2f - size_mapped=%2.f <= %.2f",
                    parts,
                    GT_MEGABYTES(suftabsize),
                    GT_MEGABYTES(largest),
                    GT_MEGABYTES(estimatedspace),
                    GT_MEGABYTES(size_mapped),
                    GT_MEGABYTES(maximumspace));
        gt_suftabparts_rc_delete(suftabparts_rc);
        return (int) parts;
      }
    }
    gt_suftabparts_rc_delete(suftabparts_rc);
  }
  gt_error_set(err,"cannot compute enhanced suffix array in at most " GT_WU
               " bytes", maximumspace);
  return -1;
}
