/*
  Copyright (c) 2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
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

#ifndef MD5SET_H
#define MD5SET_H

typedef struct GtMd5set GtMd5set;

/* a set which represents set of sequences, and allows to check
 * if a sequence or its reverse complement belongs to the set;
 * only 128bit md5 hashes are stored, not the sequence themselves
 * (a low probability of failure due to collisions exists)
 *
 * parameter nof_elements: use 0 of a lower bound if unknown
 * this parameter eliminates/reduces need to reallocate md5 table
 * increasing efficiency
 * */
GtMd5set *gt_md5set_new(unsigned long nof_elements);

void gt_md5set_delete(GtMd5set *md5set);

/* calculate the md5 hash of an upper case copy of seq
 * and add it to the set if not present*
 * if double_strand is true, the md5 hash of the reverse
 * complement is calculated too; the direct md5 of set is added
 * only if both direct and reverse complement md5 are not present
 *
 * return value:
 *
 * < 0:  error
 * 0:    md5 not found, added to set
 * 1:    md5 of seq (upcase) found, nothing added
 * 2:    md5 of rev.compl. of seq (upcase) found, nothing added
 * */
int gt_md5set_add_sequence(GtMd5set *md5set, const char* seq,
    unsigned long seqlen, bool double_strand, GtError *err);

#endif
