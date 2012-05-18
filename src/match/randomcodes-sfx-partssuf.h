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

#ifndef RANDOMCODES_SFX_PARTSSUF_H
#define RANDOMCODES_SFX_PARTSSUF_H

#include "core/logger.h"
#include "core/codetype.h"
#include "core/error_api.h"
#include "sfx-maprange.h"
#include "randomcodes-tab.h"
#include "bcktab.h"

typedef struct GtSuftabparts_rc GtSuftabparts_rc;

void gt_suftabparts_rc_showallrecords(const GtSuftabparts_rc *suftabparts_rc,
                                   bool withminmaxindex);

GtSuftabparts_rc *gt_suftabparts_rc_new(unsigned int numofparts,
                                  const GtBcktab *bcktab,
                                  const GtRandomcodestab *fct,
                                  const GtSfxmappedrangelist *sfxmrlist,
                                  unsigned long numofsuffixestoinsert,
                                  unsigned long fullspecials,
                                  GtLogger *logger);

void gt_suftabparts_rc_delete(GtSuftabparts_rc *suftabparts_rc);

GtCodetype gt_suftabparts_rc_minindex(unsigned int part,
                                   const GtSuftabparts_rc *suftabparts_rc);

GtCodetype gt_suftabparts_rc_maxindex(unsigned int part,
                                   const GtSuftabparts_rc *suftabparts_rc);

GtCodetype gt_suftabparts_rc_maxindex_last(
    const GtSuftabparts_rc *suftabparts_rc);

unsigned long gt_suftabparts_rc_offset(unsigned int part,
                                    const GtSuftabparts_rc *suftabparts_rc);

unsigned long gt_suftabparts_rc_sumofwidth(unsigned int part,
                                        const GtSuftabparts_rc *suftabparts_rc);

unsigned long gt_suftabparts_rc_widthofpart(unsigned int part,
    const GtSuftabparts_rc *suftabparts_rc);

unsigned long gt_suftabparts_rc_largest_width(
    const GtSuftabparts_rc *suftabparts_rc);

unsigned int gt_suftabparts_rc_numofparts(
    const GtSuftabparts_rc *suftabparts_rc);

unsigned long gt_suftabparts_rc_largestsizemappedpartwise(
                                       const GtSuftabparts_rc *suftabparts_rc);

double gt_suftabparts_rc_variance(const GtSuftabparts_rc *suftabparts_rc);

int gt_suftabparts_rc_fit_memlimit(size_t estimatedspace,
                                unsigned long maximumspace,
                                const GtBcktab *bcktab,
                                const GtRandomcodestab *fct,
                                const GtSfxmappedrangelist *sfxmrlist,
                                unsigned long totallength,
                                unsigned int bitsforseqnumrelpos,
                                unsigned long specialcharacters,
                                unsigned long numofsuffixestosort,
                                bool suftabuint,
                                GtError *err);

#endif
