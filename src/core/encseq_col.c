/*
  Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#include "core/class_alloc_lock.h"
#include "core/cstr_api.h"
#include "core/encseq_col.h"
#include "core/grep_api.h"
#include "core/encseq.h"
#include "core/ma.h"
#include "core/md5_seqid.h"
#include "core/seq_col_rep.h"
#include "core/seq_info_cache.h"
#include "core/undef_api.h"

struct GtEncseqCol {
  GtSeqCol parent_instance;
  GtEncseq *encseq;
  GtMD5Tab *md5_tab;
  GtSeqInfoCache *grep_cache;
};

const GtSeqColClass* gt_encseq_col_class(void);
#define gt_encseq_col_cast(SC)\
        gt_seq_col_cast(gt_encseq_col_class(), SC)

static void gt_encseq_col_delete(GtSeqCol *sc)
{
  GtEncseqCol *esc;
  esc = gt_encseq_col_cast(sc);
  if (!esc) return;
  gt_seq_info_cache_delete(esc->grep_cache);
  gt_md5_tab_delete(esc->md5_tab);
  gt_encseq_delete(esc->encseq);
}

static int gt_encseq_col_do_grep_desc(GtEncseqCol *esc, GtUword *filenum,
                                      GtUword *seqnum, GtStr *seqid,
                                      GtError *err)
{
  GtUword j;
  const GtSeqInfo *seq_info_ptr;
  GtSeqInfo seq_info;
  bool match = false;
  int had_err = 0;
  gt_error_check(err);

  gt_assert(esc->encseq && gt_encseq_has_description_support(esc->encseq));

  gt_assert(esc && filenum && seqnum && seqid);
  /* create cache */
  if (!esc->grep_cache)
    esc->grep_cache = gt_seq_info_cache_new();
  /* try to read from cache */
  seq_info_ptr = gt_seq_info_cache_get(esc->grep_cache, gt_str_get(seqid));
  if (seq_info_ptr) {
    *filenum = seq_info_ptr->filenum;
    *seqnum = seq_info_ptr->seqnum;
    return 0;
  }
  for (j = 0; !had_err && j < gt_encseq_num_of_sequences(esc->encseq); j++) {
    const char *desc;
    GtUword desc_len;
    desc = gt_encseq_description(esc->encseq, &desc_len, j);
    gt_assert(desc);
    had_err = gt_grep_nt(&match, gt_str_get(seqid), desc, desc_len, err);
    if (!had_err && match) {
      *filenum = seq_info.filenum =
                       gt_encseq_filenum(esc->encseq,
                                         gt_encseq_seqstartpos(esc->encseq, j));
      *seqnum = seq_info.seqnum =
                      j - gt_encseq_filenum_first_seqnum(esc->encseq, *filenum);
      gt_seq_info_cache_add(esc->grep_cache, gt_str_get(seqid), &seq_info);
      break;
    }
  }
  if (!had_err && !match) {
    gt_error_set(err, "no description matched sequence ID '%s'",
                 gt_str_get(seqid));
    had_err = -1;
  }
  return had_err;
}

static int gt_encseq_col_grep_desc(GtSeqCol *sc, char **seq,
                                   GtUword start, GtUword end,
                                   GtStr *seqid, GtError *err)
{
  GtUword filenum = 0, seqnum = 0;
  int had_err;
  GtEncseqCol *esc;
  esc = gt_encseq_col_cast(sc);
  gt_error_check(err);
  gt_assert(esc && seq && seqid);
  had_err = gt_encseq_col_do_grep_desc(esc, &filenum, &seqnum, seqid, err);
  if (!had_err) {
    *seq = gt_seq_col_get_sequence(sc, filenum, seqnum, start, end);
  }
  return had_err;
}

static int gt_encseq_col_grep_desc_md5(GtSeqCol *sc, const char **md5,
                                       GtStr *seqid, GtError *err)
{
  GtUword filenum = 0, seqnum = 0;
  int had_err;
  GtEncseqCol *esc;
  esc = gt_encseq_col_cast(sc);
  gt_error_check(err);
  gt_assert(esc && md5 && seqid);
  had_err = gt_encseq_col_do_grep_desc(esc, &filenum, &seqnum, seqid, err);
  if (!had_err)
    *md5 = gt_seq_col_get_md5_fingerprint(sc, seqnum, filenum);
  return had_err;
}

static int gt_encseq_col_grep_desc_sequence_length(GtSeqCol *sc,
                                                   GtUword *length,
                                                   GtStr *seqid,
                                                   GtError *err)
{
  GtUword filenum = 0, seqnum = 0;
  int had_err;
  GtEncseqCol *esc;
  esc = gt_encseq_col_cast(sc);
  gt_error_check(err);
  gt_assert(esc && length && seqid);
  had_err = gt_encseq_col_do_grep_desc(esc, &filenum, &seqnum, seqid, err);
  if (!had_err)
    *length = gt_seq_col_get_sequence_length(sc, seqnum, filenum);
  return had_err;
}

static int gt_encseq_col_md5_to_seq(GtSeqCol *sc, char **seq,
                                    GtUword start, GtUword end,
                                    GtStr *md5_seqid, GtError *err)
{
  GtUword seqnum = GT_UNDEF_UWORD;
  char seqid[GT_MD5_SEQID_HASH_LEN + 1];
  int had_err = 0;
  GtEncseqCol *esc;
  esc = gt_encseq_col_cast(sc);
  gt_error_check(err);
  gt_assert(esc && seq && start <= end && md5_seqid && err);
  gt_assert(gt_md5_seqid_has_prefix(gt_str_get(md5_seqid)));
  if (gt_str_length(md5_seqid) >= GT_MD5_SEQID_TOTAL_LEN) {
    const char *cstrseqid = gt_str_get(md5_seqid);
    if (cstrseqid[GT_MD5_SEQID_TOTAL_LEN-1] != GT_MD5_SEQID_SEPARATOR) {
      gt_error_set(err, "MD5 sequence id %s not terminated with '%c'",
                   gt_str_get(md5_seqid), GT_MD5_SEQID_SEPARATOR);
      had_err = -1;
    }
    if (!had_err) {
      strncpy(seqid, cstrseqid + GT_MD5_SEQID_PREFIX_LEN,
              GT_MD5_SEQID_HASH_LEN);
      seqid[GT_MD5_SEQID_HASH_LEN] = '\0';
    }
  }
  seqnum = gt_md5_tab_map(esc->md5_tab, seqid);
  if (seqnum != GT_UNDEF_UWORD) {
    GtUword startpos = gt_encseq_seqstartpos(esc->encseq, seqnum),
                  GT_UNUSED seqlength = gt_encseq_seqlength(esc->encseq,
                                                            seqnum);
    *seq = gt_calloc(end - start + 1, sizeof (char));
    gt_encseq_extract_decoded(esc->encseq, (char*) *seq, startpos + start,
                              startpos + end);
  } else {
    gt_error_set(err, "sequence %s not found", gt_str_get(md5_seqid));
    had_err = -1;
  }
  return had_err;
}

static int gt_encseq_col_md5_to_description(GtSeqCol *sc, GtStr *desc,
                                            GtStr *md5_seqid, GtError *err)
{
  GtUword seqnum = GT_UNDEF_UWORD;
  char seqid[GT_MD5_SEQID_HASH_LEN + 1];
  int had_err = 0;
  GtEncseqCol *esc;
  esc = gt_encseq_col_cast(sc);
  gt_error_check(err);
  gt_assert(esc && desc && md5_seqid && err);
  gt_assert(gt_md5_seqid_has_prefix(gt_str_get(md5_seqid)));
  if (gt_str_length(md5_seqid) >= GT_MD5_SEQID_TOTAL_LEN) {
    const char *cstrseqid = gt_str_get(md5_seqid);
    if (cstrseqid[GT_MD5_SEQID_TOTAL_LEN-1] != GT_MD5_SEQID_SEPARATOR) {
      gt_error_set(err, "MD5 sequence id %s not terminated with '%c'",
                   gt_str_get(md5_seqid), GT_MD5_SEQID_SEPARATOR);
      had_err = -1;
    }
    if (!had_err) {
      strncpy(seqid, cstrseqid + GT_MD5_SEQID_PREFIX_LEN,
              GT_MD5_SEQID_HASH_LEN);
      seqid[GT_MD5_SEQID_HASH_LEN] = '\0';
    }
  }
  seqnum = gt_md5_tab_map(esc->md5_tab, seqid);
  if (seqnum != GT_UNDEF_UWORD) {
    const char *cdesc;
    GtUword desc_len;
    gt_assert(seqnum < gt_encseq_num_of_sequences(esc->encseq));
    cdesc = gt_encseq_description(esc->encseq, &desc_len, seqnum);
    gt_str_append_cstr_nt(desc, cdesc, desc_len);
  } else {
    gt_error_set(err, "sequence %s not found", gt_str_get(md5_seqid));
    had_err = -1;
  }
  return had_err;
}

int gt_encseq_col_md5_to_sequence_length(GtSeqCol *sc, GtUword *len,
                                         GtStr *md5_seqid, GtError *err)
{
  GtUword seqnum = GT_UNDEF_UWORD;
  int had_err = 0;
  GtEncseqCol *esc;
  esc = gt_encseq_col_cast(sc);
  gt_error_check(err);
  gt_assert(esc && len && md5_seqid && err);
  gt_assert(gt_md5_seqid_has_prefix(gt_str_get(md5_seqid)));
  seqnum = gt_md5_tab_map(esc->md5_tab, gt_str_get(md5_seqid) +
                                          GT_MD5_SEQID_PREFIX_LEN);
  if (seqnum != GT_UNDEF_UWORD) {
    gt_assert(seqnum < gt_encseq_num_of_sequences(esc->encseq));
    *len = gt_encseq_seqlength(esc->encseq, seqnum);
  } else {
    gt_error_set(err, "sequence %s not found", gt_str_get(md5_seqid));
    had_err = -1;
  }
  return had_err;
}

static GtUword gt_encseq_col_num_of_files(const GtSeqCol *sc)
{
  const GtEncseqCol *esc;
  esc = gt_encseq_col_cast(sc);
  gt_assert(esc);
  return gt_encseq_num_of_files(esc->encseq);
}

static GtUword gt_encseq_col_num_of_seqs(const GtSeqCol *sc,
                                               GtUword filenum)
{
  GtEncseqCol *esc;
  /* XXX cache function evaluated values */
  esc = gt_encseq_col_cast(sc);
  gt_assert(esc && filenum < gt_encseq_num_of_files(esc->encseq));
  if (gt_encseq_num_of_files(esc->encseq) == 1 && filenum == 0)
    return gt_encseq_num_of_sequences(esc->encseq);
  else if (filenum == gt_encseq_num_of_files(esc->encseq) - 1) {
    return (gt_encseq_num_of_sequences(esc->encseq)
              - gt_encseq_filenum_first_seqnum(esc->encseq, filenum));
  } else {
    GtUword firstpos, nextpos;
    gt_assert(filenum < gt_encseq_num_of_files(esc->encseq) - 1);
    firstpos = gt_encseq_filenum_first_seqnum(esc->encseq, filenum);
    nextpos = gt_encseq_filenum_first_seqnum(esc->encseq, filenum + 1);
    return nextpos - firstpos;
  }
}

static const char* gt_encseq_col_get_md5_fingerprint(const GtSeqCol *sc,
                                                     GtUword filenum,
                                                     GtUword seqnum)
{
  GtEncseqCol *esc;
  esc = gt_encseq_col_cast(sc);
  gt_assert(esc && filenum < gt_encseq_num_of_files(esc->encseq));
  return gt_md5_tab_get(esc->md5_tab,
                        gt_encseq_filenum_first_seqnum(esc->encseq,
                                                       filenum) + seqnum);
}

static char* gt_encseq_col_get_sequence(const GtSeqCol *sc,
                                        GtUword filenum,
                                        GtUword seqnum,
                                        GtUword start,
                                        GtUword end)
{
  GtEncseqCol *esc;
  char *out;
  GtUword encseq_seqnum, startpos;
  esc = gt_encseq_col_cast(sc);
  gt_assert(esc && filenum < gt_encseq_num_of_files(esc->encseq));
  encseq_seqnum = gt_encseq_filenum_first_seqnum(esc->encseq, filenum) + seqnum;
  gt_assert(encseq_seqnum < gt_encseq_num_of_sequences(esc->encseq));
  gt_assert(start <= end);
  startpos = gt_encseq_seqstartpos(esc->encseq, encseq_seqnum);
  out = gt_calloc(end - start + 1, sizeof (char));
  gt_encseq_extract_decoded(esc->encseq, out, startpos + start, startpos + end);
  return out;
}

static char* gt_encseq_col_get_description(const GtSeqCol *sc,
                                           GtUword filenum,
                                           GtUword seqnum)
{
  GtEncseqCol *esc;
  const char *desc;
  GtUword encseq_seqnum, desclen;
  esc = gt_encseq_col_cast(sc);
  gt_assert(esc && filenum < gt_encseq_num_of_files(esc->encseq));
  encseq_seqnum = gt_encseq_filenum_first_seqnum(esc->encseq, filenum) + seqnum;
  gt_assert(encseq_seqnum < gt_encseq_num_of_sequences(esc->encseq));
  desc = gt_encseq_description(esc->encseq, &desclen, encseq_seqnum);
  gt_assert(desc && desclen > 0);
  return gt_cstr_dup_nt(desc, desclen);;
}

static GtUword gt_encseq_col_get_sequence_length(const GtSeqCol *sc,
                                                       GtUword filenum,
                                                       GtUword seqnum)
{
  GtEncseqCol *esc;
  GtUword encseq_seqnum;
  esc = gt_encseq_col_cast(sc);
  gt_assert(esc && filenum < gt_encseq_num_of_files(esc->encseq));
  encseq_seqnum = gt_encseq_filenum_first_seqnum(esc->encseq, filenum) + seqnum;
  return gt_encseq_seqlength(esc->encseq, encseq_seqnum);
}

const GtSeqColClass* gt_encseq_col_class(void)
{
  static const GtSeqColClass *esc_class = NULL;
  gt_class_alloc_lock_enter();
  if (!esc_class) {
      esc_class = gt_seq_col_class_new(sizeof (GtEncseqCol),
                                       gt_encseq_col_delete,
                                       gt_encseq_col_grep_desc,
                                       gt_encseq_col_grep_desc_md5,
                                       gt_encseq_col_grep_desc_sequence_length,
                                       gt_encseq_col_md5_to_seq,
                                       gt_encseq_col_md5_to_description,
                                       gt_encseq_col_md5_to_sequence_length,
                                       gt_encseq_col_num_of_files,
                                       gt_encseq_col_num_of_seqs,
                                       gt_encseq_col_get_md5_fingerprint,
                                       gt_encseq_col_get_sequence,
                                       gt_encseq_col_get_description,
                                       gt_encseq_col_get_sequence_length);
  }
  gt_class_alloc_lock_leave();
  return esc_class;
}

GtSeqCol* gt_encseq_col_new(GtEncseq *encseq, GtError *err)
{
  GtSeqCol *sc;
  GtEncseqCol *esc;
  gt_error_check(err);
  gt_assert(encseq);
  if (!gt_encseq_has_md5_support(encseq)) {
    gt_error_set(err, "encoded sequence has no MD5 support");
    return NULL;
  }
  sc = gt_seq_col_create(gt_encseq_col_class());
  esc = gt_encseq_col_cast(sc);
  esc->md5_tab = gt_encseq_get_md5_tab(encseq, err);
  gt_assert(esc->md5_tab);
  esc->encseq = gt_encseq_ref(encseq);
  return sc;
}
