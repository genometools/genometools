/*
  Copyright (c) 2007-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef REGION_MAPPING_H
#define REGION_MAPPING_H

#include "core/str.h"

/* maps a sequence-region to a sequence file */
typedef struct GtRegionMapping GtRegionMapping;

GtRegionMapping* gt_region_mapping_new_mapping(GtStr *mapping_filename,
                                               GtError*);
GtRegionMapping* gt_region_mapping_new_seqfile(GtStr *sequence_filename,
                                               bool usedesc);
GtRegionMapping* gt_region_mapping_ref(GtRegionMapping*);
int              gt_region_mapping_get_raw_sequence(GtRegionMapping*,
                                                    const char **raw,
                                                    GtStr *seqid, GtError*);
int              gt_region_mapping_get_raw_sequence_length(GtRegionMapping*,
                                                           unsigned long*,
                                                           GtStr *seqid,
                                                           GtError*);
void             gt_region_mapping_delete(GtRegionMapping*);

#endif
