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

#define FROMCODE2SPECIALCODE(CODE,NUMOFCHARS)\
                            (((NUMOFCHARS) == 4U)\
                            ? ((CODE) >> 2)\
                            : (((CODE) - ((NUMOFCHARS)-1)) / (NUMOFCHARS)))

typedef enum
{
  GtSfxGtBitsequence,
  GtSfxunsignedlong,
  GtSfxuint32_t
} GtSfxmappedrangetype;

typedef struct GtSfxmappedrange GtSfxmappedrange;

typedef struct
{
  unsigned long mapoffset, mapend;
} GtMappedrange;

void gt_mapped_lbrange_get(GtMappedrange *range,
                           size_t sizeofbasetype,
                           unsigned long pagesize,
                           unsigned long mincode,
                           unsigned long maxcode);

void gt_mapped_csrange_get(GtMappedrange *range,
                           unsigned int padoffset,
                           size_t sizeofbasetype,
                           unsigned long numofallcodes,
                           unsigned int numofchars,
                           unsigned long pagesize,
                           GtCodetype mincode,
                           GtCodetype maxcode);

unsigned int gt_Sfxmappedrange_padoffset(size_t sizeofbasetype,
                                         unsigned long offset);

GtSfxmappedrange *gt_Sfxmappedrange_new(void **usedptrptr,
                                        unsigned long numofindexes,
                                        GtSfxmappedrangetype type,
                                        const char *tablename,
                                        GtLogger *logger,
                                        GtError *err);

void gt_Sfxmappedrange_make_writable(GtSfxmappedrange *sfxmappedrange);

void *gt_Sfxmappedrange_map(GtSfxmappedrange *sfxmappedrange,
                            unsigned int part,
                            unsigned long minindex,
                            unsigned long maxindex,
                            GtLogger *logger);

int gt_Sfxmappedrange_delete(GtSfxmappedrange *sfxmappedrange,
                             GtLogger *logger,GtError *err);

int gt_unlink_possibly_with_error(const char *filename,GtLogger *logger,
                                  GtError *err);

#endif
