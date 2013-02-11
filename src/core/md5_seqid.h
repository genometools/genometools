/*
  Copyright (c) 2010      Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c)      2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c)      2013 Center for Bioinformatics, University of Hamburg

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

#ifndef MD5_SEQID_H
#define MD5_SEQID_H

#include <stdbool.h>
#include "core/error_api.h"

#define GT_MD5_SEQID_PREFIX      "md5:"
#define GT_MD5_SEQID_PREFIX_LEN  4
#define GT_MD5_SEQID_HASH_LEN    32
#define GT_MD5_SEQID_TOTAL_LEN   (GT_MD5_SEQID_PREFIX_LEN + \
                                  GT_MD5_SEQID_HASH_LEN + 1)
#define GT_MD5_SEQID_SEPARATOR   ':'

/* Returns <true> if <seqid> has the prefix used to denote MD5 sequence IDs,
   <false> otherwise. */
bool gt_md5_seqid_has_prefix(const char *seqid);

/* Compares \0-terminated seqid strings <id_a> and <id_b> (similarly to
   strcmp(3)), but is aware of MD5 prefixes. That is, if both seqids have MD5
   prefixes, only the MD5 prefixes will be compared. If at least one seqid has
   no MD5 prefix, the sequd without the prefix will sort before the other
   one. */
int  gt_md5_seqid_cmp_seqids(const char *id_a, const char *id_b);

int  gt_md5_seqid_unit_test(GtError *err);

#endif
