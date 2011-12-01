/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#ifndef SFX_MAPRANGE_H
#define SFX_MAPRANGE_H

#include <stdbool.h>
#include <stdlib.h>
#include "core/str_api.h"
#include "core/codetype.h"

typedef enum
{
  GtSfxGtBitsequence,
  GtSfxunsignedlong,
  GtSfxuint32_t
} GtSfxmappedrangetype;

typedef struct GtSfxmappedrange GtSfxmappedrange;

void *gt_Sfxmappedrange_map_entire(GtSfxmappedrange *sfxmappedrange,
                                   GtError *err);

typedef void (*GtSfxmappedrangetransformfunc)(unsigned long *,
                                              unsigned long *,
                                              const void *);

GtSfxmappedrange *gt_Sfxmappedrange_new(const char *tablename,
                                        unsigned long numofentries,
                                        GtSfxmappedrangetype type,
                                        GtSfxmappedrangetransformfunc
                                          transformfunc,
                                        const void *transformfunc_data);

void gt_Sfxmappedrange_storetmp_ulong(GtSfxmappedrange *sfxmappedrange,
                                      unsigned long **usedptrptr,
                                      bool writable);

void gt_Sfxmappedrange_storetmp_uint32(GtSfxmappedrange *sfxmappedrange,
                                       uint32_t **usedptrptr,
                                       bool writable);

void gt_Sfxmappedrange_storetmp_bitsequence(GtSfxmappedrange *sfxmappedrange,
                                            GtBitsequence **usedptrptr,
                                            bool writable);

void gt_Sfxmappedrange_usetmp(GtSfxmappedrange *sfxmappedrange,
                              const GtStr *tmpfilename,
                              void **usedptrptr,
                              unsigned long numofentries,
                              bool writable);

void *gt_Sfxmappedrange_map(GtSfxmappedrange *sfxmappedrange,
                            unsigned long minindex,
                            unsigned long maxindex);

void gt_Sfxmappedrange_unmap(GtSfxmappedrange *sfxmappedrange);

size_t gt_Sfxmappedrange_size_entire(const GtSfxmappedrange *sfxmappedrange);

unsigned long gt_Sfxmappedrange_size_mapped(const GtSfxmappedrange
                                                      *sfxmappedrange,
                                            unsigned long minindex,
                                            unsigned long maxindex);

void gt_Sfxmappedrange_checkindex(const GtSfxmappedrange *sfxmappedrange,
                                  GT_UNUSED unsigned long idx);

void gt_Sfxmappedrange_delete(GtSfxmappedrange *sfxmappedrange);

typedef struct GtSfxmappedrangelist GtSfxmappedrangelist;

GtSfxmappedrangelist *gt_Sfxmappedrangelist_new(void);

void gt_Sfxmappedrangelist_add(GtSfxmappedrangelist *sfxmrlist,
                               GtSfxmappedrange *sfxmappedrange);

unsigned long gt_Sfxmappedrangelist_size_mapped(
                                         const GtSfxmappedrangelist *sfxmrlist,
                                         unsigned long minindex,
                                         unsigned long maxindex);

void gt_Sfxmappedrangelist_delete(GtSfxmappedrangelist *sfxmrlist);

unsigned long gt_Sfxmappedrangelist_size_entire(
                                       const GtSfxmappedrangelist *sfxmrlist);

#endif
