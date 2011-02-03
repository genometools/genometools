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

#include "core/fa.h"
#include "core/fileutils_api.h"
#include "core/hashmap_api.h"
#include "core/ma.h"
#include "core/md5_fingerprint.h"
#include "core/md5tab.h"
#include "core/undef.h"
#include "core/xansi_api.h"

struct GtMD5Tab{
  FILE *fingerprints_file; /* used to lock the memory mapped fingerprints */
  char *fingerprints; /* holds memory mapped fingerprints */
  char **md5_fingerprints;
  unsigned long num_of_md5s;
  bool owns_md5s,
       has_map;
  GtHashmap *md5map; /* maps md5 to index */
};

static bool read_fingerprints(GtMD5Tab *md5tab,
                              GtStr *fingerprints_filename)
{
  bool reading_succeeded = true;
  size_t len;
  gt_assert(md5tab && fingerprints_filename);
  /* open file */
  if (gt_file_exists(gt_str_get(fingerprints_filename))) {
    md5tab->fingerprints_file = gt_fa_xfopen(gt_str_get(fingerprints_filename),
                                             "r");
    gt_fa_lock_shared(md5tab->fingerprints_file);
    md5tab->fingerprints = gt_fa_xmmap_read(gt_str_get(fingerprints_filename),
                                            &len);
    if (len != md5tab->num_of_md5s * 33) {
      gt_fa_xmunmap(md5tab->fingerprints);
      md5tab->fingerprints = NULL;
      gt_fa_unlock(md5tab->fingerprints_file);
      gt_fa_xfclose(md5tab->fingerprints_file);
      md5tab->fingerprints_file = NULL;
      reading_succeeded = false;
    }
  }
  else
    reading_succeeded = false;
  /* reading file (each line contains a single MD5 sum) */
  if (reading_succeeded) {
    unsigned long i;
    for (i = 0; i < md5tab->num_of_md5s; i++)
      md5tab->md5_fingerprints[i] = md5tab->fingerprints + i * 33;
  }
  return reading_succeeded;
}

static void add_fingerprints(char **md5_fingerprints, const void *seqs,
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
                               GtStr *fingerprints_filename)
{
  FILE *fingerprints_file;
  gt_assert(md5_fingerprints && num_of_md5s && fingerprints_filename);
  fingerprints_file = gt_fa_xfopen(gt_str_get(fingerprints_filename), "w");
  gt_fa_lock_exclusive(fingerprints_file);
  dump_md5_fingerprints(md5_fingerprints, num_of_md5s, fingerprints_file);
  gt_fa_unlock(fingerprints_file);
  gt_fa_xfclose(fingerprints_file);
}

GtMD5Tab* gt_md5tab_new(const char *sequence_file, const void *seqs,
                        GtGetSeqFunc get_seq, GtGetSeqLenFunc get_seq_len,
                        unsigned long num_of_seqs, bool use_cache_file,
                        bool build_map)
{
  GtMD5Tab *md5tab;
  bool reading_succeeded = false;
  GtStr *fingerprints_filename;
  unsigned long i;
  gt_assert(sequence_file && seqs && get_seq && get_seq_len);
  md5tab = gt_calloc(1, sizeof *md5tab);
  md5tab->md5_fingerprints = gt_calloc(num_of_seqs, sizeof (char*));
  md5tab->num_of_md5s = num_of_seqs;
  fingerprints_filename = gt_str_new_cstr(sequence_file);
  gt_str_append_cstr(fingerprints_filename, GT_MD5TAB_FILE_SUFFIX);
  if (use_cache_file && gt_file_exists(gt_str_get(fingerprints_filename)) &&
      !gt_file_is_newer(sequence_file, gt_str_get(fingerprints_filename))) {
    /* only try to read the fingerprint file if the sequence file was not
       modified in the meantime */
    reading_succeeded = read_fingerprints(md5tab, fingerprints_filename);
  }
  if (!reading_succeeded) {
    add_fingerprints(md5tab->md5_fingerprints, seqs, get_seq, get_seq_len,
                     num_of_seqs);
    md5tab->owns_md5s = true;
    if (use_cache_file) {
      write_fingerprints(md5tab->md5_fingerprints, md5tab->num_of_md5s,
                         fingerprints_filename);
    }
  }
  gt_str_delete(fingerprints_filename);
  /* fill md5 map, if necessary */
  if (build_map) {
    md5tab->has_map = true;
    md5tab->md5map = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
    for (i = 0; i < md5tab->num_of_md5s; i++) {
      gt_hashmap_add(md5tab->md5map, (void*) md5tab->md5_fingerprints[i],
                     (void*) (i + 1));
    }
  }
  return md5tab;
}

void gt_md5tab_delete(GtMD5Tab *md5tab)
{
  unsigned long i;
  if (!md5tab) return;
  gt_fa_xmunmap(md5tab->fingerprints);
  gt_fa_unlock(md5tab->fingerprints_file);
  gt_fa_xfclose(md5tab->fingerprints_file);
  gt_hashmap_delete(md5tab->md5map);
  if (md5tab->owns_md5s) {
    for (i = 0; i < md5tab->num_of_md5s; i++)
      gt_free(md5tab->md5_fingerprints[i]);
  }
  gt_free(md5tab->md5_fingerprints);
  gt_free(md5tab);
}

const char* gt_md5tab_get(const GtMD5Tab *md5tab, unsigned long idx)
{
  gt_assert(md5tab && idx < md5tab->num_of_md5s);
  return md5tab->md5_fingerprints[idx];
}

unsigned long gt_md5tab_map(const GtMD5Tab *md5tab, const char *md5)
{
  const char *value;
  gt_assert(md5tab && md5);
  gt_assert(md5tab->has_map);
  value = gt_hashmap_get(md5tab->md5map, md5);
  if (value)
    return ((unsigned long) value) - 1;
  return GT_UNDEF_ULONG;
}
