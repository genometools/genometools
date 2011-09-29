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

#ifndef SFX_PARTSSUF_H
#define SFX_PARTSSUF_H

#include "core/logger.h"
#include "core/codetype.h"
#include "core/error_api.h"
#include "sfx-maprange.h"
#include "firstcodes-tab.h"
#include "bcktab.h"

typedef struct GtSuftabparts GtSuftabparts;

void gt_suftabparts_showallrecords(const GtSuftabparts *suftabparts,
                                   bool withminmaxindex);

GtSuftabparts *gt_suftabparts_new(unsigned int numofparts,
                                  const GtBcktab *bcktab,
                                  const GtFirstcodestab *fct,
                                  const GtSfxmappedrangelist *sfxmrlist,
                                  unsigned long numofsuffixestoinsert,
                                  unsigned long fullspecials,
                                  GtLogger *logger);

void gt_suftabparts_delete(GtSuftabparts *suftabparts);

GtCodetype gt_suftabparts_minindex(unsigned int part,
                                   const GtSuftabparts *suftabparts);

GtCodetype gt_suftabparts_maxindex(unsigned int part,
                                   const GtSuftabparts *suftabparts);

GtCodetype gt_suftabparts_maxindex_last(const GtSuftabparts *suftabparts);

unsigned long gt_suftabparts_offset(unsigned int part,
                                    const GtSuftabparts *suftabparts);

unsigned long gt_suftabparts_sumofwidth(unsigned int part,
                                        const GtSuftabparts *suftabparts);

unsigned long gt_suftabparts_widthofpart(unsigned int part,
                                         const GtSuftabparts *suftabparts);

unsigned long gt_suftabparts_largest_width(const GtSuftabparts *suftabparts);

unsigned int gt_suftabparts_numofparts(const GtSuftabparts *suftabparts);

unsigned long gt_suftabparts_largestsizemappedpartwise(
                                       const GtSuftabparts *suftabparts);

double gt_suftabparts_variance(const GtSuftabparts *suftabparts);

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
                                GtError *err);

#endif
