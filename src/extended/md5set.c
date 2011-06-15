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

#include <string.h>
#include <ctype.h>
#include "md5.h"
#include "core/ma.h"
#include "extended/reverse.h"
#include "core/hashtable.h"
#include "core/safearith.h"
#include "extended/md5set.h"

#define GT_MD5SET_MD5SIZE (sizeof (unsigned char) * 16)

struct GtMd5set
{
  GtHashtable *ht;
  HashElemInfo info;
  char* buffer;
  unsigned long bufsize;
};

htsize_t gt_md5set_md5_hash(const void *md5)
{
  htsize_t hash;

  gt_assert(sizeof (htsize_t) < GT_MD5SET_MD5SIZE);
  memcpy(&hash, md5, sizeof (htsize_t));
  return hash;
}

int gt_md5set_md5_cmp(const void *md5A, const void *md5B)
{
  return memcmp(md5A, md5B, GT_MD5SET_MD5SIZE);
}

GtMd5set *gt_md5set_new(void)
{
  GtMd5set *md5set;
  md5set = gt_malloc(sizeof (GtMd5set));
  md5set->info.keyhash = gt_md5set_md5_hash;
  md5set->info.free_op.free_elem = NULL;
  md5set->info.elem_size = 16;
  md5set->info.cmp = gt_md5set_md5_cmp;
  md5set->info.table_data = NULL;
  md5set->info.table_data_free = NULL;
  md5set->ht = gt_hashtable_new(md5set->info);
  md5set->buffer = NULL;
  md5set->bufsize = 0;
  return md5set;
}

void gt_md5set_delete(GtMd5set *md5set)
{
  if (md5set != NULL)
  {
    gt_hashtable_delete(md5set->ht);
    gt_free(md5set->buffer);
    gt_free(md5set);
  }
}

static void gt_md5set_prepare_buffer(GtMd5set *md5set, unsigned long bufsize)
{
  gt_assert(md5set != NULL);

  if (md5set->buffer == NULL)
  {
    md5set->buffer = gt_malloc(sizeof (char) * bufsize);
    md5set->bufsize = bufsize;
  }
  else if (md5set->bufsize < bufsize)
  {
    md5set->buffer = gt_realloc(md5set->buffer, sizeof (char) * bufsize);
    md5set->bufsize = bufsize;
  }
}

#define GT_MD5SET_NOT_FOUND  0
#define GT_MD5SET_FOUND      1
#define GT_MD5SET_RC_FOUND   2

int gt_md5set_add_sequence(GtMd5set *md5set, const char* seq,
    unsigned long seqlen, bool double_strand, GtError *err)
{
  char md5sum[16], md5sum_rc[16];
  unsigned long i;
  int retval = 0;

  gt_assert(md5set != NULL);
  gt_assert(md5set->ht != NULL);

  gt_md5set_prepare_buffer(md5set, seqlen);
  for (i = 0; i < seqlen; i++)
    md5set->buffer[i] = toupper(seq[i]);

  md5(md5set->buffer, gt_safe_cast2long(seqlen), md5sum);

  if (gt_hashtable_get(md5set->ht, md5sum))
  {
    return GT_MD5SET_FOUND;
  }
  else
  {
    gt_hashtable_add(md5set->ht, md5sum);
    if (double_strand)
    {
      retval = gt_reverse_complement(md5set->buffer, seqlen, err);
      if (retval != 0)
      {
        gt_assert(retval < 0);
        return retval;
      }
      md5(md5set->buffer, gt_safe_cast2long(seqlen), md5sum_rc);
      if (gt_hashtable_get(md5set->ht, md5sum_rc))
      {
        return GT_MD5SET_RC_FOUND;
      }
    }
    return GT_MD5SET_NOT_FOUND;
  }
}
