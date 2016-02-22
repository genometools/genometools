/*
  Copyright (c) 2011      Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c)      2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2011      Center for Bioinformatics, University of Hamburg
  Copyright (c) 2010-2012 Center for Bioinformatics, University of Hamburg

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

#ifndef SAM_ALIGNMENT_REP_H
#define SAM_ALIGNMENT_REP_H

/* The contents of this file is to be considered private implementation detail.
*/

#include <samtools/sam.h>
#include "core/alphabet_api.h"

struct GtSamAlignment{
  bam1_t       *s_alignment;
  GtAlphabet   *alphabet;
  GtUchar      *seq_buffer,
               *qual_buffer;
  GtUword s_bufsize,
                q_bufsize,
                rightmost;
};

#endif
