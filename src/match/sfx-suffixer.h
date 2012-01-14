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

#ifndef SFX_SUFFIXER_H
#define SFX_SUFFIXER_H

#include "core/encseq_api.h"
#include "core/error.h"
#include "core/logger.h"
#include "core/readmode.h"
#include "core/timer_api.h"
#include "core/intbits.h"
#include "core/codetype.h"
#include "extended/uint64hashtable.h"
#include "firstcodes-buf.h"
#include "sfx-suffixgetset.h"
#include "sfx-strategy.h"

typedef struct Sfxiterator Sfxiterator;

int gt_Sfxiterator_delete(Sfxiterator *sfi,GtError *err);

Sfxiterator *gt_Sfxiterator_new(const GtEncseq *encseq,
                                GtReadmode readmode,
                                unsigned int prefixlength,
                                unsigned int numofparts,
                                unsigned long maximumspace,
                                const Sfxstrategy *sfxstrategy,
                                GtTimer *sfxprogress,
                                bool withprogressbar,
                                GtLogger *logger,
                                GtError *err);

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
                                GtError *err);

const GtSuffixsortspace *gt_Sfxiterator_next(unsigned long *numberofsuffixes,
                                             bool *specialsuffixes,
                                             Sfxiterator *sfi);

int gt_Sfxiterator_postsortfromstream(Sfxiterator *sfi,
                                      const GtStr *indexname,
                                      GtError *err);

int gt_Sfxiterator_bcktab2file(FILE *fp,Sfxiterator *sfi,GtError *err);

unsigned long gt_Sfxiterator_longest(const Sfxiterator *sfi);

GtCodetype gt_kmercode_at_firstpos(const GtTwobitencoding *twobitencoding,
                                   unsigned int kmersize);

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
                                   void *processkmerspecialinfo);

#endif
