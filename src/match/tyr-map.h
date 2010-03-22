/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#ifndef TYR_MAP_H
#define TYR_MAP_H

#include "core/str_api.h"
#include "core/error_api.h"
#include "core/symboldef.h"
#include "core/defined-types.h"

typedef struct Tyrindex Tyrindex;

Tyrindex *tyrindex_new(const GtStr *tyrindexname,GtError *err);
const GtUchar *tyrindex_mertable(const Tyrindex *tyrindex);
const GtUchar *tyrindex_lastmer(const Tyrindex *tyrindex);
unsigned long tyrindex_merbytes(const Tyrindex *tyrindex);
unsigned int tyrindex_alphasize(const Tyrindex *tyrindex);
unsigned long tyrindex_mersize(const Tyrindex *tyrindex);
bool tyrindex_isempty(const Tyrindex *tyrindex);
void tyrindex_show(const Tyrindex *tyrindex);
void tyrindex_delete(Tyrindex **tyrindexptr);
/*@null@*/ const GtUchar *tyrindex_binmersearch(const Tyrindex *tyrindex,
                                              unsigned long offset,
                                              const GtUchar *key,
                                              const GtUchar *leftbound,
                                              const GtUchar *rightbound);
void tyrindex_check(const Tyrindex *tyrindex);
int determinetyrbckpfxlen(unsigned int *prefixlength,
                          const Tyrindex *tyrindex,
                          const Definedunsignedint *callprefixlength,
                          GtError *err);
unsigned long tyrindex_ptr2number(const Tyrindex *tyrindex,
                                  const GtUchar *result);
typedef struct Tyrcountinfo Tyrcountinfo;

Tyrcountinfo *tyrcountinfo_new(const Tyrindex *tyrindex,
                               const GtStr *tyrindexname,
                               GtError *err);
unsigned long tyrcountinfo_get(const Tyrcountinfo *tyrcountinfo,
                               unsigned long mernumber);
void tyrcountinfo_delete(Tyrcountinfo **tyrcountinfoptr);

#endif
