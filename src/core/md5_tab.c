/*
  Copyright (c) 2006-2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/fa.h"
#include "core/fileutils_api.h"
#include "core/hashmap_api.h"
#include "core/ma.h"
#include "core/md5_fingerprint_api.h"
#include "core/md5_tab.h"
#include "core/undef_api.h"
#include "core/xansi_api.h"

struct GtMD5Tab{
  FILE *fingerprints_file; /* used to lock the memory mapped fingerprints */
  char *fingerprints; /* holds memory mapped fingerprints */
  char **md5_fingerprints;
  unsigned long num_of_md5s,
                reference_count;
  bool owns_md5s;
  GtHashmap *md5map; /* maps md5 to index */
};

static bool read_fingerprints(GtMD5Tab *md5_tab,
                              const char *fingerprints_filename,
                              bool use_file_locking)
{
  bool reading_succeeded = true;
  size_t len;
  gt_assert(md5_tab && fingerprints_filename);
  /* open file */
  gt_assert(gt_file_exists(fingerprints_filename));
  if (use_file_locking) {
    md5_tab->fingerprints_file = gt_fa_xfopen(fingerprints_filename, "r");
    gt_fa_lock_shared(md5_tab->fingerprints_file);
  }
  md5_tab->fingerprints = gt_fa_xmmap_read(fingerprints_filename, &len);
  if (len != md5_tab->num_of_md5s * 33) {
    gt_fa_xmunmap(md5_tab->fingerprints);
    md5_tab->fingerprints = NULL;
    gt_fa_unlock(md5_tab->fingerprints_file);
    gt_fa_xfclose(md5_tab->fingerprints_file);
    md5_tab->fingerprints_file = NULL;
    reading_succeeded = false;
  }
  return reading_succeeded;
}

static void add_fingerprints(char **md5_fingerprints, void *seqs,
                             GtGetSeqFunc get_seq, GtGetSeqLenFunc get_seq_len,
                             unsigned long num_of_seqs)
{
  unsigned long i;
  gt_assert(md5_fingerprints && seqs && get_seq && get_seq_len);
  for (i = 0; i < num_of_seqs; i++) {
    md5_fingerprints[i] = gt_md5_fingerprint(get_seq(seqs, i),
                                             get_seq_len(seqs, i));
  }
}

static void dump_md5_fingerprints(char **md5_fingerprints,
                                  unsigned long num_of_md5s, FILE *outfp)
{
  unsigned long i;
  gt_assert(md5_fingerprints && num_of_md5s && outfp);
  for (i = 0; i < num_of_md5s; i++) {
    gt_xfputs(md5_fingerprints[i], outfp);
    gt_xfputc('\0', outfp);
  }
}

static void write_fingerprints(char **md5_fingerprints,
                               unsigned long num_of_md5s,
                               GtStr *fingerprints_filename,
                               bool use_file_locking)
{
  FILE *fingerprints_file;
  gt_assert(md5_fingerprints && num_of_md5s && fingerprints_filename);
  fingerprints_file = gt_fa_xfopen(gt_str_get(fingerprints_filename), "w");
  if (use_file_locking)
    gt_fa_lock_exclusive(fingerprints_file);
  dump_md5_fingerprints(md5_fingerprints, num_of_md5s, fingerprints_file);
  if (use_file_locking)
    gt_fa_unlock(fingerprints_file);
  gt_fa_xfclose(fingerprints_file);
}

GtMD5Tab* gt_md5_tab_new(const char *sequence_file, void *seqs,
                         GtGetSeqFunc get_seq, GtGetSeqLenFunc get_seq_len,
                         unsigned long num_of_seqs, bool use_cache_file,
                         bool use_file_locking)
{
  GtMD5Tab *md5_tab;
  bool reading_succeeded = false;
  GtStr *fingerprints_filename;
  gt_assert(sequence_file && seqs && get_seq && get_seq_len);
  md5_tab = gt_calloc(1, sizeof *md5_tab);
  md5_tab->num_of_md5s = num_of_seqs;
  fingerprints_filename = gt_str_new_cstr(sequence_file);
  gt_str_append_cstr(fingerprints_filename, GT_MD5_TAB_FILE_SUFFIX);
  if (use_cache_file && gt_file_exists(gt_str_get(fingerprints_filename)) &&
      !gt_file_is_newer(sequence_file, gt_str_get(fingerprints_filename))) {
    /* only try to read the fingerprint file if the sequence file was not
       modified in the meantime */
    reading_succeeded = read_fingerprints(md5_tab,
                                          gt_str_get(fingerprints_filename),
                                          use_file_locking);
  }
  if (!reading_succeeded) {
    md5_tab->md5_fingerprints = gt_calloc(num_of_seqs, sizeof (char*));
    add_fingerprints(md5_tab->md5_fingerprints, seqs, get_seq, get_seq_len,
                     num_of_seqs);
    md5_tab->owns_md5s = true;
    if (use_cache_file) {
      write_fingerprints(md5_tab->md5_fingerprints, md5_tab->num_of_md5s,
                         fingerprints_filename, use_file_locking);
    }
  }
  gt_str_delete(fingerprints_filename);
  return md5_tab;
}

GtMD5Tab* gt_md5_tab_new_from_cache_file(const char *cache_file,
                                         unsigned long num_of_seqs,
                                         bool use_file_locking,
                                         GtError *err)
{
  GtMD5Tab *md5_tab;
  bool reading_succeeded = false;
  gt_assert(cache_file);
  gt_error_check(err);

  md5_tab = gt_calloc(1, sizeof *md5_tab);
  md5_tab->num_of_md5s = num_of_seqs;
  if (gt_file_exists(cache_file)) {
    reading_succeeded = read_fingerprints(md5_tab,
                                          cache_file,
                                          use_file_locking);
  }
  if (!reading_succeeded) {
    gt_free(md5_tab);
    gt_error_set(err, "could not read fingerprints file \"%s\" or "
                      "invalid file contents", cache_file);
    return NULL;
  }
  md5_tab->owns_md5s = false;
  return md5_tab;
}

GtMD5Tab* gt_md5_tab_ref(GtMD5Tab *md5_tab)
{
  if (!md5_tab) return NULL;
  md5_tab->reference_count++;
  return md5_tab;
}

void gt_md5_tab_delete(GtMD5Tab *md5_tab)
{
  unsigned long i;
  if (!md5_tab) return;
  if (md5_tab->reference_count) {
    md5_tab->reference_count--;
    return;
  }
  gt_fa_xmunmap(md5_tab->fingerprints);
  gt_fa_unlock(md5_tab->fingerprints_file);
  gt_fa_xfclose(md5_tab->fingerprints_file);
  gt_hashmap_delete(md5_tab->md5map);
  if (md5_tab->owns_md5s) {
    for (i = 0; i < md5_tab->num_of_md5s; i++)
      gt_free(md5_tab->md5_fingerprints[i]);
    gt_free(md5_tab->md5_fingerprints);
  }
  gt_free(md5_tab);
}

const char* gt_md5_tab_get(const GtMD5Tab *md5_tab, unsigned long idx)
{
  gt_assert(md5_tab && idx < md5_tab->num_of_md5s);
  if (md5_tab->owns_md5s)
    return md5_tab->md5_fingerprints[idx];
 return md5_tab->fingerprints + idx * 33;
}

static void build_md5map(GtMD5Tab *md5_tab)
{
  unsigned long i;
  gt_assert(md5_tab);
  md5_tab->md5map = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  for (i = 0; i < md5_tab->num_of_md5s; i++) {
    gt_hashmap_add(md5_tab->md5map, (void*) gt_md5_tab_get(md5_tab, i),
                   (void*) (i + 1));
  }
}

unsigned long gt_md5_tab_map(GtMD5Tab *md5_tab, const char *md5)
{
  const char *value;
  gt_assert(md5_tab && md5);
  if (!md5_tab->md5map)
    build_md5map(md5_tab);
  gt_assert(md5_tab->md5map);
  value = gt_hashmap_get(md5_tab->md5map, md5);
  if (value)
    return ((unsigned long) value) - 1;
  return GT_UNDEF_ULONG;
}

unsigned long gt_md5_tab_size(const GtMD5Tab *md5_tab)
{
  gt_assert(md5_tab);
  return md5_tab->num_of_md5s;
}
