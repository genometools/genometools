/*
  Copyright (c) 2010      Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c)      2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c)      2012 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "core/assert_api.h"
#include "core/ensure.h"
#include "core/md5_seqid.h"
#include "core/undef_api.h"

bool gt_md5_seqid_has_prefix(const char *seqid)
{
  gt_assert(seqid);
  return !strncmp(seqid, GT_MD5_SEQID_PREFIX, GT_MD5_SEQID_PREFIX_LEN);
}

int gt_md5_seqid_cmp_seqids(const char *id_a, const char *id_b)
{
  int md5_status = 0;
  gt_assert(id_a && id_b);

  if (id_a == id_b)
    return 0;
  if (gt_md5_seqid_has_prefix(id_a))
    md5_status += 1;
  if (gt_md5_seqid_has_prefix(id_b))
    md5_status += 2;

  switch (md5_status) {
    case 0:
      return strcmp(id_a, id_b);
    case 1:
      return 1;
    case 2:
      return -1;
    case 3:
      return strncmp(id_a, id_b, GT_MD5_SEQID_TOTAL_LEN);
    default:
      gt_assert(false);   /* cannot happen */
  }
  return 0; /* cannot happen */
}

int gt_md5_seqid_unit_test(GtError *err) {
  int had_err = 0;
  const char *seqid1 = "foo",
             *seqid1_diffptr = "_foo",
             *seqid1_md5 = "md5:d3b07384d113edec49eaa6238ad5ff00:foo",
             *seqid1_wrongmd5 = "md5:c157a79031e1c40f85931829bc5fc552:foo";
  gt_error_check(err);

  gt_ensure(had_err, !gt_md5_seqid_has_prefix(seqid1));
  gt_ensure(had_err, gt_md5_seqid_has_prefix(seqid1_md5));
  gt_ensure(had_err, gt_md5_seqid_cmp_seqids(seqid1, seqid1) == 0);
  gt_ensure(had_err, gt_md5_seqid_cmp_seqids(seqid1, seqid1_diffptr+1) == 0);
  gt_ensure(had_err, gt_md5_seqid_cmp_seqids(seqid1, seqid1_md5) == -1);
  gt_ensure(had_err, gt_md5_seqid_cmp_seqids(seqid1_md5, seqid1) == 1);
  gt_ensure(had_err, gt_md5_seqid_cmp_seqids(seqid1_md5, seqid1_md5) == 0);
  gt_ensure(had_err, gt_md5_seqid_cmp_seqids(seqid1_md5, seqid1_wrongmd5) > 0);

  return had_err;
}
