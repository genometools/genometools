/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2013 Ole Eigenbrod <ole.eigenbrod@gmx.de>
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

#ifndef SFX_MAPPEDSTR_H
#define SFX_MAPPEDSTR_H

#include "core/error_api.h"
#include "core/str_array_api.h"
#include "core/encseq_api.h"
#include "core/codetype.h"
#include "intcode-def.h"

typedef struct
{
  bool definedspecialposition;
  unsigned int specialposition;
  GtCodetype code;
} GtKmercode;

typedef struct GtKmercodeiterator GtKmercodeiterator;

/*@notnull@*/
GtKmercodeiterator *gt_kmercodeiterator_encseq_new(const GtEncseq *encseq,
                                                   GtReadmode readmode,
                                                   unsigned int kmersize,
                                                   GtUword startpos);

GtUword             gt_kmercodeiterator_encseq_get_currentpos(
                                                GtKmercodeiterator *kmercodeit);

void                gt_kmercodeiterator_encseq_set_currentpos(
                                                 GtKmercodeiterator *kmercodeit,
                                                 GtUword position);

void                gt_kmercodeiterator_reset(
                                         GtKmercodeiterator *kmercodeiterator,
                                         GtReadmode readmode, GtUword startpos);

bool                gt_kmercodeiterator_encseq_isspecial(
                                                GtKmercodeiterator *kmercodeit);
bool                gt_kmercodeiterator_encseq_isexhausted(
                                                GtKmercodeiterator *kmercodeit);
void                gt_kmercodeiterator_encseq_setexhausted(
                                                 GtKmercodeiterator *kmercodeit,
                                                 bool exhausted);

const GtKmercode*   gt_kmercodeiterator_encseq_next(
                                          GtKmercodeiterator *kmercodeiterator);

const GtKmercode*   gt_kmercodeiterator_encseq_nonspecial_next(
                                          GtKmercodeiterator *kmercodeiterator);

GtKmercodeiterator* gt_kmercodeiterator_filetab_new(
                                                  const GtStrArray *filenametab,
                                                  unsigned int numofchars,
                                                  unsigned int kmersize,
                                                  const GtUchar *symbolmap,
                                                  bool plainformat,
                                                  GtError *err);

int                 gt_kmercodeiterator_filetab_next(
                                           const GtKmercode **kmercodeptr,
                                           GtKmercodeiterator *kmercodeiterator,
                                           GtError *err);

bool                gt_kmercodeiterator_inputexhausted(
                                    const GtKmercodeiterator *kmercodeiterator);

void                gt_kmercodeiterator_delete(
                                          GtKmercodeiterator *kmercodeiterator);

void                getencseqkmers(const GtEncseq *encseq,
                                   GtReadmode readmode,
                                   unsigned int kmersize,
                                   void(*processkmercode)(void *,
                                                          GtUword,
                                                          const GtKmercode *),
                                   void *processkmercodeinfo);

#endif
