/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef REGIONMAPPING_H
#define REGIONMAPPING_H

#include "libgtcore/str.h"

/* maps a sequence-region to a sequence file */
typedef struct RegionMapping RegionMapping;

RegionMapping* regionmapping_new_mapping(Str *mapping_filename, Env*);
RegionMapping* regionmapping_new_seqfile(Str *sequence_filename, Env*);
int            regionmapping_get_raw_sequence(RegionMapping*, const char**,
                                              Str *seqid, Env*);
int            regionmapping_get_raw_sequence_length(RegionMapping*,
                                                     unsigned long*,
                                                     Str *seqid, Env*);
void           regionmapping_delete(RegionMapping*, Env*);

#endif
