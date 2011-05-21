/*
  Copyright (c) 2008-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008      Center for Bioinformatics, University of Hamburg

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

#ifndef REGION_FACTORY_H
#define REGION_FACTORY_H

#include "extended/node_visitor_api.h"
#include "gth/input.h"

typedef struct GthRegionFactory GthRegionFactory;

GthRegionFactory* gth_region_factory_new(bool use_desc_ranges);
void              gth_region_factory_delete(GthRegionFactory*);
/* Use <region_factory> to produce sequence regions for each genomic sequence in
   <input> and save them in the array <nodes>. */
void              gth_region_factory_save(GthRegionFactory *region_factory,
                                          GtArray *nodes, GthInput *input);
GtStr*            gth_region_factory_get_seqid(GthRegionFactory*,
                                               unsigned long filenum,
                                               unsigned long seqnum);
long              gth_region_factory_offset(GthRegionFactory*,
                                            unsigned long filenum,
                                            unsigned long seqnum);

#endif
