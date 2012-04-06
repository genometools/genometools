/*
  Copyright (c) 2010, 2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#include "core/bioseq.h"
#include "core/bioseq_collection.h"
#include "core/cstr_api.h"
#include "core/grep_api.h"
#include "core/hashmap_api.h"
#include "core/ma.h"
#include "core/md5_seqid.h"
#include "core/undef_api.h"

struct GtBioseqCollection {
  GtBioseq **bioseqs;
  unsigned long num_of_seqfiles;
  GtHashmap *grep_cache;
};

GtBioseqCollection* gt_bioseq_collection_new(GtStrArray *sequence_files,
                                             GtError *err)
{
  GtBioseqCollection *bsc;
  unsigned long i;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(sequence_files);
  gt_assert(gt_str_array_size(sequence_files));
  bsc = gt_calloc(1, sizeof *bsc);
  bsc->num_of_seqfiles = gt_str_array_size(sequence_files);
  bsc->bioseqs = gt_calloc(bsc->num_of_seqfiles, sizeof (GtBioseq*));
  for (i = 0; !had_err && i < bsc->num_of_seqfiles; i++) {
    bsc->bioseqs[i] = gt_bioseq_new(gt_str_array_get(sequence_files, i), err);
    if (!bsc->bioseqs[i])
      had_err = -1;
  }
  if (had_err) {
    gt_bioseq_collection_delete(bsc);
    return NULL;
  }
  return bsc;
}

void gt_bioseq_collection_delete(GtBioseqCollection *bsc)
{
  unsigned long i;
  if (!bsc) return;
  gt_hashmap_delete(bsc->grep_cache);
  for (i = 0; i < bsc->num_of_seqfiles; i++)
    gt_bioseq_delete(bsc->bioseqs[i]);
  gt_free(bsc->bioseqs);
  gt_free(bsc);
}

typedef struct {
  unsigned long filenum,
                seqnum;
} SeqInfo;

static int grep_desc(GtBioseqCollection *bsc, unsigned long *filenum,
                     unsigned long *seqnum, GtStr *seqid, GtError *err)
{
  unsigned long i, j;
  SeqInfo *seq_info;
  bool match;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(bsc && filenum && seqnum && seqid);
  /* create cache */
  if (!bsc->grep_cache) {
    bsc->grep_cache = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                     gt_free_func);
  }
  /* try to read from cache */
  if ((seq_info = gt_hashmap_get(bsc->grep_cache, gt_str_get(seqid)))) {
    *filenum = seq_info->filenum;
    *seqnum = seq_info->seqnum;
    return 0;
  }
  for (i = 0; !had_err && i < bsc->num_of_seqfiles; i++) {
    GtBioseq *bioseq = bsc->bioseqs[i];
    for (j = 0; !had_err && j < gt_bioseq_number_of_sequences(bioseq); j++) {
      const char *desc = gt_bioseq_get_description(bioseq, j);
      had_err = gt_grep(&match, gt_str_get(seqid), desc, err);
      if (!had_err && match) {
        *filenum = i;
        *seqnum = j;
        /* cache results */
        seq_info = gt_malloc(sizeof *seq_info);
        seq_info->filenum = i;
        seq_info->seqnum = j;
        gt_hashmap_add(bsc->grep_cache, gt_cstr_dup(gt_str_get(seqid)),
                       seq_info);
        break;
      }
    }
    if (match)
      break;
  }
  if (!had_err && !match) {
    gt_error_set(err, "no description matched sequence ID '%s'",
                 gt_str_get(seqid));
    had_err = -1;
  }
  return had_err;
}

int gt_bioseq_collection_grep_desc(GtBioseqCollection *bsc, const char **rawseq,
                                   unsigned long *length, GtStr *seqid,
                                   GtError *err)
{
  unsigned long filenum = 0, seqnum = 0;
  int had_err;
  gt_error_check(err);
  gt_assert(bsc && rawseq && length && seqid);
  had_err = grep_desc(bsc, &filenum, &seqnum, seqid, err);
  if (!had_err) {
    *rawseq = gt_bioseq_get_sequence(bsc->bioseqs[filenum], seqnum);
    *length = gt_bioseq_get_sequence_length(bsc->bioseqs[filenum], seqnum);
  }
  return had_err;
}

int gt_bioseq_collection_grep_desc_md5(GtBioseqCollection *bsc,
                                       const char **md5, GtStr *seqid,
                                       GtError *err)
{
  unsigned long filenum = 0, seqnum = 0;
  int had_err;
  gt_error_check(err);
  gt_assert(bsc && md5 && seqid);
  had_err = grep_desc(bsc, &filenum, &seqnum, seqid, err);
  if (!had_err)
    *md5 = gt_bioseq_get_md5_fingerprint(bsc->bioseqs[filenum], seqnum);
  return had_err;
}

int gt_bioseq_collection_md5_to_seq(GtBioseqCollection *bsc, const char **seq,
                                    unsigned long *length, GtStr *md5_seqid,
                                    GtError *err)
{
  unsigned long i, seqnum = GT_UNDEF_ULONG;
  int had_err = 0;
  GtBioseq *bioseq;
  gt_error_check(err);
  gt_assert(bsc && seq && length && md5_seqid && err);
  gt_assert(gt_md5_seqid_has_prefix(gt_str_get(md5_seqid)));
  /* XXX: extract method */
  for (i = 0; i < bsc->num_of_seqfiles; i++) {
    bioseq = bsc->bioseqs[i];
    seqnum = gt_bioseq_md5_to_index(bioseq, gt_str_get(md5_seqid) +
                                    GT_MD5_SEQID_PREFIX_LEN);
    if (seqnum != GT_UNDEF_ULONG)
      break;
  }
  if (seqnum != GT_UNDEF_ULONG) {
    *seq = gt_bioseq_get_sequence(bioseq, seqnum);
    *length = gt_bioseq_get_sequence_length(bioseq, seqnum);
  }
  else  {
    gt_error_set(err, "sequence %s not found", gt_str_get(md5_seqid));
    had_err = -1;
  }
  return had_err;
}
int gt_bioseq_collection_md5_to_description(GtBioseqCollection *bsc,
                                            GtStr *desc, GtStr *md5_seqid,
                                            GtError *err)
{
  unsigned long i, seqnum = GT_UNDEF_ULONG;
  int had_err = 0;
  GtBioseq *bioseq;
  gt_error_check(err);
  gt_assert(bsc && desc && md5_seqid && err);
  gt_assert(gt_md5_seqid_has_prefix(gt_str_get(md5_seqid)));
  /* XXX: extract method */
  for (i = 0; i < bsc->num_of_seqfiles; i++) {
    bioseq = bsc->bioseqs[i];
    seqnum = gt_bioseq_md5_to_index(bioseq, gt_str_get(md5_seqid) +
                                    GT_MD5_SEQID_PREFIX_LEN);
    if (seqnum != GT_UNDEF_ULONG)
      break;
  }
  if (seqnum != GT_UNDEF_ULONG)
    gt_str_append_cstr(desc, gt_bioseq_get_description(bioseq, seqnum));
  else  {
    gt_error_set(err, "sequence %s not found", gt_str_get(md5_seqid));
    had_err = -1;
  }
  return had_err;
}
