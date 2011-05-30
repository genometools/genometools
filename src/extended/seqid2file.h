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

#ifndef SEQID2FILE_H
#define SEQID2FILE_H

#include "core/option_api.h"
#include "extended/region_mapping_api.h"

typedef struct GtSeqid2FileInfo GtSeqid2FileInfo;

GtSeqid2FileInfo* gt_seqid2file_info_new(void);
void              gt_seqid2file_info_delete(GtSeqid2FileInfo*);

/* add the options -seqfile and -regionmapping to the given option parser */
void              gt_seqid2file_register_options(GtOptionParser*,
                                                 GtSeqid2FileInfo*);

GtRegionMapping*  gt_seqid2file_region_mapping_new(GtSeqid2FileInfo*,
                                                   GtError*);

#endif
