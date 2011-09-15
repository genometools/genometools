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
#include "core/logger_api.h"
#include "core/codetype.h"

typedef enum
{
  GtSfxGtBitsequence,
  GtSfxunsignedlong,
  GtSfxuint32_t
} GtSfxmappedrangetype;

typedef struct GtSfxmappedrange GtSfxmappedrange;
typedef struct GtSfxmappedrangelist GtSfxmappedrangelist;

GtSfxmappedrangelist *gt_Sfxmappedrangelist_new(size_t size_of_elem);

void gt_Sfxmappedrangelist_add(GtSfxmappedrangelist *sfxmrlist,
                               GtSfxmappedrange *sfxmappedrange);

unsigned long gt_Sfxmappedrange_size_of_part(
                                         const GtSfxmappedrangelist *sfxmrlist,
                                         unsigned long minindex,
                                         unsigned long maxindex);

void gt_Sfxmappedrangelist_delete(GtSfxmappedrangelist *sfxmrlist);

void *gt_Sfxmappedrange_map_entire(GtSfxmappedrange *sfxmappedrange,
                                   GtError *err);

size_t gt_Sfxmappedrange_size_entire(const GtSfxmappedrange *sfxmappedrange);

GtSfxmappedrange *gt_Sfxmappedrange_new(unsigned long numofentries,
                                        GtSfxmappedrangetype type,
                                        unsigned long(*transformfunc)(
                                            unsigned long,unsigned int),
                                        unsigned int transformfunc_data);

int gt_Sfxmappedrange_enhance(GtSfxmappedrange *sfxmappedrange,
                              void **usedptrptr,
                              bool writable,
                              const char *tablename,
                              GtLogger *logger,
                              GtError *err);

void *gt_Sfxmappedrange_map(GtSfxmappedrange *sfxmappedrange,
                            unsigned int part,
                            unsigned long minindex,
                            unsigned long maxindex,
                            GtLogger *logger);

unsigned long gt_Sfxmappedrange_size_mapped(const GtSfxmappedrange
                                              *sfxmappedrange,
                                            unsigned long minindex,
                                            unsigned long maxindex);

void gt_Sfxmappedrange_delete(GtSfxmappedrange *sfxmappedrange,
                              GtLogger *logger);

#endif
