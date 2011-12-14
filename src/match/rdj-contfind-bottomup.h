/*
  Copyright (c) 2010-2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2010-2011 Center for Bioinformatics, University of Hamburg

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

#ifndef RDJ_CONTFIND_BOTTOMUP_H
#define RDJ_CONTFIND_BOTTOMUP_H

#include "match/esa-seqread.h"
#include "core/intbits.h"

/*
* Bottom-up dfs-based esa algorithm to identify contained reads.
*
* A read r is contained in r' iff:
*   - r is a proper substring of r' or revcompl(r'):
*   - r == r' or r == revcompl(r') and r'.seqnum < r.seqnum
*
* Input:
* - GtBitsequence *contained:
*   - size must be 1 bit/read
*   - should be initialized with 0
* - firstrevcompl:
*   - no mirror, single-strand mode     -> 0
*   - virtual mirror                    -> nofreads
*   - real mirror                       -> nofseqs / 2
* - read_length:
*   - fixed length       -> read length  [modulo operations are used]
*   - variable length    -> 0            [a bit table is created]
*
* Output:
* - return value: number of contained reads
* - GtBitsequence *contained:
*   bits are set for each read respecting the definition
*/
unsigned long gt_contfind_bottomup(Sequentialsuffixarrayreader *ssar,
    bool show_progressbar, GtBitsequence *contained,
    unsigned long firstrevcompl, unsigned long read_length);

#endif
