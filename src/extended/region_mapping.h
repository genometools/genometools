/*
  Copyright (c) 2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

#include "extended/region_mapping_api.h"

/* Like <gt_region_mapping_new_encseq()>, but interpreting a seqid of "seq<n>"
   as sequence <n> in <encseq>. Used as an access method for legacy LTRharvest
   output files. */
GtRegionMapping* gt_region_mapping_new_encseq_seqno(GtEncseq *encseq);
/* Enables matching only at the beginning of sequence descriptions up to the
   first whitespace */
void             gt_region_mapping_enable_match_desc_start(GtRegionMapping *rm);

#endif
