/*
  Copyright (c) 2006-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef MD5_TAB_H
#define MD5_TAB_H

#define GT_MD5_TAB_FILE_SUFFIX ".md5"

typedef struct GtMD5Tab GtMD5Tab;

typedef const char*   (*GtGetSeqFunc)(const void *seqs, unsigned long index);
typedef unsigned long (*GtGetSeqLenFunc)(const void *seqs, unsigned long index);

/* Create a new MD5 table object for sequences contained in <sequence_file>.
   The sequences have to be stored in <seqs> (<num_of_seqs> many) and have to be
   accessible via the functions <get_seq> and <get_seq_len>. If <use_cache_file>
   is <true>, the MD5 sums are read from a cache file (named
   "<sequence_file><GT_MD5TAB_FILE_SUFFIX>"), if it exists or written to it, if
   it doesn't exist. If <use_cache_file> is <false>, no cache file is read or
   written. */
GtMD5Tab*     gt_md5_tab_new(const char *sequence_file, const void *seqs,
                             GtGetSeqFunc get_seq, GtGetSeqLenFunc get_seq_len,
                             unsigned long num_of_seqs, bool use_cache_file);
void          gt_md5_tab_delete(GtMD5Tab *md5_tab);
/* Return the MD5 sum for sequence <index>. */
const char*   gt_md5_tab_get(const GtMD5Tab*, unsigned long index);
/* Map <md5> back to sequence index. */
unsigned long gt_md5_tab_map(GtMD5Tab*, const char *md5);
unsigned long gt_md5_tab_size(const GtMD5Tab*);

#endif
