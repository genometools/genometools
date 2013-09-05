/*
  Copyright (c) 2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#ifndef RDJ_CONTIG_INFO_H
#define RDJ_CONTIG_INFO_H

#include <stdint.h>

typedef struct {
  uint64_t deg :16;
  uint64_t ptr :48;
} GtContigEdgesLink;

typedef struct {
  GtUword  depth;
  GtUword  length;
  float          firstread_copynum;
  float          internal_copynum;
  float          lastread_copynum;
  float          astat;
} GtContigDepthInfo;

typedef struct {
  GtUword junction_num;
  uint32_t contig_num;
  unsigned int length     :31;
  unsigned int firstnode  :1;
} GtContigJunctionInfo;

#endif
