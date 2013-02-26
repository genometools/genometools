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
#include "core/bioseq_col.h"
#include "core/class_alloc_lock.h"
#include "core/cstr_api.h"
#include "core/grep_api.h"
#include "core/hashmap_api.h"
#include "core/ma.h"
#include "core/md5_seqid.h"
#include "core/seq_col_rep.h"
#include "core/seq_info_cache.h"
#include "core/undef_api.h"

struct GtBioseqCol {
  GtSeqCol parent_instance;
  GtBioseq **bioseqs;
  unsigned long num_of_seqfiles;
  GtSeqInfoCache *grep_cache;
};

const GtSeqColClass* gt_bioseq_col_class(void);
#define gt_bioseq_col_cast(SC)\
        gt_seq_col_cast(gt_bioseq_col_class(), SC)

static void gt_bioseq_col_delete(GtSeqCol *sc)
{
  unsigned long i;
  GtBioseqCol *bsc;
  bsc = gt_bioseq_col_cast(sc);
  if (!bsc) return;
  gt_seq_info_cache_delete(bsc->grep_cache);
  for (i = 0; i < bsc->num_of_seqfiles; i++)
    gt_bioseq_delete(bsc->bioseqs[i]);
  gt_free(bsc->bioseqs);
}

static int grep_desc(GtBioseqCol *bsc, unsigned long *filenum,
                     unsigned long *seqnum, GtStr *seqid, GtError *err)
{
  unsigned long i, j;
  const GtSeqInfo *seq_info_ptr;
  GtSeqInfo seq_info;
  bool match;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(bsc && filenum && seqnum && seqid);
  /* create cache */
  if (!bsc->grep_cache)
    bsc->grep_cache = gt_seq_info_cache_new();
  /* try to read from cache */
  seq_info_ptr = gt_seq_info_cache_get(bsc->grep_cache, gt_str_get(seqid));
  if (seq_info_ptr) {
    *filenum = seq_info_ptr->filenum;
    *seqnum = seq_info_ptr->seqnum;
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
        seq_info.filenum = i;
        seq_info.seqnum = j;
        gt_seq_info_cache_add(bsc->grep_cache, gt_str_get(seqid), &seq_info);
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

static int gt_bioseq_col_grep_desc(GtSeqCol *sc, char **seq,
                                   unsigned long start, unsigned long end,
                                   GtStr *seqid, GtError *err)
{
  unsigned long filenum = 0, seqnum = 0, seqlength;
  int had_err;
  GtBioseqCol *bsc;
  bsc = gt_bioseq_col_cast(sc);
  gt_error_check(err);
  gt_assert(bsc && seq && seqid);
  had_err = grep_desc(bsc, &filenum, &seqnum, seqid, err);
  if (!had_err) {
    seqlength = gt_bioseq_get_sequence_length(bsc->bioseqs[filenum], seqnum);
    if (start > seqlength - 1 || end > seqlength - 1) {
      had_err = -1;
      gt_error_set(err, "trying to extract range %lu-%lu on sequence "
                         "``%s''which is not covered by that sequence (only "
                         "%lu characters in size). Has the sequence-region "
                         "to sequence mapping been defined correctly?",
                         start, end, gt_str_get(seqid), seqlength);
    }
  }
  if (!had_err) {
    *seq = gt_bioseq_get_sequence_range(bsc->bioseqs[filenum], seqnum,
                                        start, end);
  }
  return had_err;
}

static int gt_bioseq_col_grep_desc_md5(GtSeqCol *sc, const char **md5,
                                       GtStr *seqid, GtError *err)
{
  unsigned long filenum = 0, seqnum = 0;
  int had_err;
  GtBioseqCol *bsc;
  bsc = gt_bioseq_col_cast(sc);
  gt_error_check(err);
  gt_assert(bsc && md5 && seqid);
  had_err = grep_desc(bsc, &filenum, &seqnum, seqid, err);
  if (!had_err)
    *md5 = gt_bioseq_get_md5_fingerprint(bsc->bioseqs[filenum], seqnum);
  return had_err;
}

static int gt_bioseq_col_grep_desc_sequence_length(GtSeqCol *sc,
                                                   unsigned long *length,
                                                   GtStr *seqid,
                                                   GtError *err)
{
  unsigned long filenum = 0, seqnum = 0;
  int had_err;
  GtBioseqCol *bsc;
  bsc = gt_bioseq_col_cast(sc);
  gt_error_check(err);
  gt_assert(bsc && length && seqid);
  had_err = grep_desc(bsc, &filenum, &seqnum, seqid, err);
  if (!had_err)
    *length = gt_bioseq_get_sequence_length(bsc->bioseqs[filenum], seqnum);
  return had_err;
}

static int gt_bioseq_col_md5_to_seq(GtSeqCol *sc, char **seq,
                                    unsigned long start, unsigned long end,
                                    GtStr *md5_seqid, GtError *err)
{
  unsigned long i, seqnum = GT_UNDEF_ULONG;
  char *seqid = NULL;
  int had_err = 0;
  GtBioseq *bioseq = NULL;
  GtBioseqCol *bsc;
  bsc = gt_bioseq_col_cast(sc);
  gt_error_check(err);
  gt_assert(bsc && seq && md5_seqid && err);
  gt_assert(gt_md5_seqid_has_prefix(gt_str_get(md5_seqid)));
  /* performance hack to avoid string duplication */
  if (gt_str_length(md5_seqid) >= GT_MD5_SEQID_TOTAL_LEN) {
    seqid = gt_str_get(md5_seqid);
    if (seqid[GT_MD5_SEQID_TOTAL_LEN-1] != GT_MD5_SEQID_SEPARATOR) {
      gt_error_set(err, "MD5 sequence id %s not terminated with '%c'",
                   gt_str_get(md5_seqid), GT_MD5_SEQID_SEPARATOR);
      had_err = -1;
    }
    if (!had_err)
      seqid[GT_MD5_SEQID_TOTAL_LEN-1] = '\0';
  }
  for (i = 0; !had_err && i < bsc->num_of_seqfiles; i++) {
    bioseq = bsc->bioseqs[i];
    seqnum = gt_bioseq_md5_to_index(bioseq, gt_str_get(md5_seqid) +
                                    GT_MD5_SEQID_PREFIX_LEN);
    if (seqnum != GT_UNDEF_ULONG)
      break;
  }
  if (!had_err) {
    if (gt_str_length(md5_seqid) >= GT_MD5_SEQID_TOTAL_LEN)
      seqid[GT_MD5_SEQID_TOTAL_LEN-1] = GT_MD5_SEQID_SEPARATOR;
    if (seqnum != GT_UNDEF_ULONG) {
      *seq = gt_bioseq_get_sequence_range(bioseq, seqnum, start, end);
    }
    else {
      gt_error_set(err, "sequence %s not found", gt_str_get(md5_seqid));
      had_err = -1;
    }
  }
  return had_err;
}

static int gt_bioseq_col_md5_to_description(GtSeqCol *sc, GtStr *desc,
                                            GtStr *md5_seqid, GtError *err)
{
  unsigned long i, seqnum = GT_UNDEF_ULONG;
  int had_err = 0;
  GtBioseq *bioseq;
  GtBioseqCol *bsc;
  bsc = gt_bioseq_col_cast(sc);
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

int gt_bioseq_col_md5_to_sequence_length(GtSeqCol *sc, unsigned long *len,
                                         GtStr *md5_seqid, GtError *err)
{
  unsigned long i, seqnum = GT_UNDEF_ULONG;
  int had_err = 0;
  GtBioseq *bioseq;
  GtBioseqCol *bsc;
  bsc = gt_bioseq_col_cast(sc);
  gt_error_check(err);
  gt_assert(bsc && len && md5_seqid && err);
  gt_assert(gt_md5_seqid_has_prefix(gt_str_get(md5_seqid)));
  for (i = 0; i < bsc->num_of_seqfiles; i++) {
    bioseq = bsc->bioseqs[i];
    seqnum = gt_bioseq_md5_to_index(bioseq, gt_str_get(md5_seqid) +
                                    GT_MD5_SEQID_PREFIX_LEN);
    if (seqnum != GT_UNDEF_ULONG)
      break;
  }
  if (seqnum != GT_UNDEF_ULONG)
    *len = gt_bioseq_get_sequence_length(bioseq, seqnum);
  else  {
    gt_error_set(err, "sequence %s not found", gt_str_get(md5_seqid));
    had_err = -1;
  }
  return had_err;
}

static unsigned long gt_bioseq_col_num_of_files(const GtSeqCol *sc)
{
  const GtBioseqCol *bsc;
  bsc = gt_bioseq_col_cast(sc);
  gt_assert(bsc);
  return bsc->num_of_seqfiles;
}

static unsigned long gt_bioseq_col_num_of_seqs(const GtSeqCol *sc,
                                               unsigned long filenum)
{
  GtBioseqCol *bsc;
  bsc = gt_bioseq_col_cast(sc);
  gt_assert(bsc && filenum < bsc->num_of_seqfiles);
  return gt_bioseq_number_of_sequences(bsc->bioseqs[filenum]);
}

static const char* gt_bioseq_col_get_md5_fingerprint(const GtSeqCol *sc,
                                                     unsigned long filenum,
                                                     unsigned long seqnum)
{
  GtBioseqCol *bsc;
  bsc = gt_bioseq_col_cast(sc);
  gt_assert(bsc && filenum < bsc->num_of_seqfiles);
  return gt_bioseq_get_md5_fingerprint(bsc->bioseqs[filenum], seqnum);
}

static char* gt_bioseq_col_get_sequence(const GtSeqCol *sc,
                                        unsigned long filenum,
                                        unsigned long seqnum,
                                        unsigned long start,
                                        unsigned long end)
{
  GtBioseqCol *bsc;
  bsc = gt_bioseq_col_cast(sc);
  gt_assert(bsc && filenum < bsc->num_of_seqfiles);
  return gt_bioseq_get_sequence_range(bsc->bioseqs[filenum], seqnum, start,
                                      end);
}

static char* gt_bioseq_col_get_description(const GtSeqCol *sc,
                                           unsigned long filenum,
                                           unsigned long seqnum)
{
  GtBioseqCol *bsc;
  bsc = gt_bioseq_col_cast(sc);
  gt_assert(bsc && filenum < bsc->num_of_seqfiles);
  return gt_cstr_dup(gt_bioseq_get_description(bsc->bioseqs[filenum], seqnum));
}

static unsigned long gt_bioseq_col_get_sequence_length(const GtSeqCol *sc,
                                                       unsigned long filenum,
                                                       unsigned long seqnum)
{
  GtBioseqCol *bsc;
  bsc = gt_bioseq_col_cast(sc);
  gt_assert(bsc && filenum < bsc->num_of_seqfiles);
  return gt_bioseq_get_sequence_length(bsc->bioseqs[filenum], seqnum);
}

const GtSeqColClass* gt_bioseq_col_class(void)
{
  static const GtSeqColClass *bsc_class = NULL;
  gt_class_alloc_lock_enter();
  if (!bsc_class) {
      bsc_class = gt_seq_col_class_new(sizeof (GtBioseqCol),
                                       gt_bioseq_col_delete,
                                       gt_bioseq_col_grep_desc,
                                       gt_bioseq_col_grep_desc_md5,
                                       gt_bioseq_col_grep_desc_sequence_length,
                                       gt_bioseq_col_md5_to_seq,
                                       gt_bioseq_col_md5_to_description,
                                       gt_bioseq_col_md5_to_sequence_length,
                                       gt_bioseq_col_num_of_files,
                                       gt_bioseq_col_num_of_seqs,
                                       gt_bioseq_col_get_md5_fingerprint,
                                       gt_bioseq_col_get_sequence,
                                       gt_bioseq_col_get_description,
                                       gt_bioseq_col_get_sequence_length);
  }
  gt_class_alloc_lock_leave();
  return bsc_class;
}

GtSeqCol* gt_bioseq_col_new(GtStrArray *sequence_files, GtError *err)
{
  GtSeqCol *sc;
  GtBioseqCol *bsc;
  unsigned long i;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(sequence_files);
  gt_assert(gt_str_array_size(sequence_files));
  sc = gt_seq_col_create(gt_bioseq_col_class());
  bsc = gt_bioseq_col_cast(sc);
  bsc->num_of_seqfiles = gt_str_array_size(sequence_files);
  bsc->bioseqs = gt_calloc(bsc->num_of_seqfiles, sizeof (GtBioseq*));
  for (i = 0; !had_err && i < bsc->num_of_seqfiles; i++) {
    bsc->bioseqs[i] = gt_bioseq_new(gt_str_array_get(sequence_files, i), err);
    if (!bsc->bioseqs[i])
      had_err = -1;
  }
  if (had_err) {
    gt_bioseq_col_delete(sc);
    return NULL;
  }
  return sc;
}
