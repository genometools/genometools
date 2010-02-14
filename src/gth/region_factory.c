/*
  Copyright (c) 2008-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008      Center for Bioinformatics, University of Hamburg

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

#include "core/cstr.h"
#include "core/cstr_table.h"
#include "core/basename_api.h"
#include "core/ma.h"
#include "core/hashmap_api.h"
#include "core/parseutils.h"
#include "core/undef.h"
#include "core/unused_api.h"
#include "core/xansi.h"
#include "gth/region_factory.h"

typedef struct {
  unsigned long num_of_files,
                *num_of_sequences;
  GtStr ***store;
  unsigned long **offsets;
} SeqidStore;

static SeqidStore* seqid_store_new(GthInput *input)
{
  SeqidStore *ss;
  unsigned long i, j;
  gt_assert(input);
  ss = gt_malloc(sizeof *ss);
  ss->num_of_files = gth_input_num_of_gen_files(input);
  ss->num_of_sequences = gt_calloc(ss->num_of_files, sizeof (unsigned long));
  /* allocate room for store */
  ss->store = gt_calloc(ss->num_of_files, sizeof *ss->store);
  for (i = 0; i < ss->num_of_files; i++) {
    ss->num_of_sequences[i] = gth_input_num_of_gen_seqs(input, i);
    ss->store[i] = gt_calloc(ss->num_of_sequences[i], sizeof **ss->store);
  }
  /* allocate room for offsets */
  ss->offsets = gt_malloc(ss->num_of_files * sizeof *ss->offsets);
  for (i = 0; i < ss->num_of_files; i++)
    ss->offsets[i] = gt_malloc(ss->num_of_sequences[i] * sizeof **ss->offsets);
  /* initialize offsets to undefined values */
  for (i = 0; i < ss->num_of_files; i++) {
    for (j = 0; j < ss->num_of_sequences[i]; j++)
      ss->offsets[i][j] = GT_UNDEF_ULONG;
  }
  return ss;
}

static void seqid_store_delete(SeqidStore *ss)
{
  unsigned long i, j;
  if (!ss) return;
  for (i = 0; i < ss->num_of_files; i++)
    gt_free(ss->offsets[i]);
  gt_free(ss->offsets);
  for (i = 0; i < ss->num_of_files; i++) {
    for (j = 0; j < ss->num_of_sequences[i]; j++)
      gt_str_delete(ss->store[i][j]);
    gt_free(ss->store[i]);
  }
  gt_free(ss->store);
  gt_free(ss->num_of_sequences);
  gt_free(ss);
}

static void seqid_store_add(SeqidStore *ss, unsigned long filenum,
                            unsigned long seqnum, GtStr *seqid,
                            unsigned long offset)
{
  gt_assert(ss && seqid);
  gt_assert(gt_str_length(seqid)); /* is not empty */
  gt_assert(filenum < ss->num_of_files);
  gt_assert(seqnum < ss->num_of_sequences[filenum]);
  gt_assert(!ss->store[filenum][seqnum]); /* is unused */
  ss->store[filenum][seqnum] = gt_str_clone(seqid);
  ss->offsets[filenum][seqnum] = offset == GT_UNDEF_ULONG ? 1 : offset;
}

static GtStr* seqid_store_get(SeqidStore *ss, unsigned long filenum,
                              unsigned long seqnum)
{
  GtStr *seqid;
  gt_assert(ss);
  gt_assert(filenum < ss->num_of_files);
  gt_assert(seqnum < ss->num_of_sequences[filenum]);
  gt_assert(ss->store[filenum][seqnum]); /* is used */
  seqid = ss->store[filenum][seqnum];
  gt_assert(gt_str_length(seqid)); /* is not empty */
  return seqid;
}

static unsigned long seqid_store_offset(SeqidStore *ss, unsigned long filenum,
                                        unsigned long seqnum)
{
  unsigned long offset;
  gt_assert(ss);
  gt_assert(filenum < ss->num_of_files);
  gt_assert(seqnum < ss->num_of_sequences[filenum]);
  offset = ss->offsets[filenum][seqnum];
  gt_assert(offset != GT_UNDEF_ULONG); /* is defined */
  return offset;
}

struct GthRegionFactory{
  GtCstrTable *used_seqids;
  bool factory_was_used,
       use_desc_ranges;
  SeqidStore *seqid_store;
};

GthRegionFactory* gth_region_factory_new(bool use_desc_ranges)
{
  GthRegionFactory *srf = gt_calloc(1, sizeof *srf);
  srf->used_seqids = gt_cstr_table_new();
  srf->use_desc_ranges = use_desc_ranges;
  return srf;
}

void gth_region_factory_delete(GthRegionFactory *srf)
{
  if (!srf) return;
  gt_cstr_table_delete(srf->used_seqids);
  seqid_store_delete(srf->seqid_store);
  gt_free(srf);
}

/* Range descriptions have the folowing format: III:1000001..2000000
   That is, the part between ':' and '..' denotes the offset. */
static unsigned long parse_description_range(const char *description)
{
  unsigned long i, desclen, offset;
  char *desc;
  gt_assert(description);
  desc = gt_xstrdup(description);
  desclen = strlen(desc);
  /* find ':' */
  for (i = 0; i < desclen; i++) {
    if (desc[i] == ':')
      break;
  }
  if (i == desclen) {
    /* no ':' found */
    gt_free(desc);
    return GT_UNDEF_ULONG;
  }
  desc += i + 1;
  /* find '..' */
  i = 0;
  while (desc[i] != '\0') {
    if (desc[i-1] == '.' && desc[i] == '.')
      break;
    i++;
  }
  if (desc[i] == '\0') {
    /* no '..' found */
    gt_free(desc);
    return GT_UNDEF_ULONG;
  }
  /* parse range */
  gt_assert(desc[i-1] == '.' && desc[i] == '.');
  desc[i-1] = '\0';
  if (gt_parse_ulong(&offset, desc)) {
    /* parsing failed */
    gt_free(desc);
    return GT_UNDEF_ULONG;
  }
  gt_free(desc);
  return offset;
}

static void make_sequence_region(GtHashmap *sequence_regions,
                                 GtStr *sequenceid,
                                 GthRegionFactory *srf,
                                 GthInput *input,
                                 unsigned long filenum,
                                 unsigned long seqnum)
{
  GtGenomeNode *sr = NULL;
  GtRange range;
  unsigned long offset = GT_UNDEF_ULONG;
  gt_assert(sequence_regions && sequenceid && srf && input);
  if (gth_input_use_substring_spec(input)) {
    range.start = gth_input_genomic_substring_from(input);
    range.end   = gth_input_genomic_substring_to(input);
  }
  else {
    range = gth_input_get_relative_genomic_range(input, filenum, seqnum);
  }
  if (srf->use_desc_ranges) {
    GtStr *description = gt_str_new();
    gth_input_get_genomic_description(input, description, filenum, seqnum);
    offset = parse_description_range(gt_str_get(description));
    gt_str_delete(description);
  }
  if (offset != GT_UNDEF_ULONG)
    range = gt_range_offset(&range, offset);
  else
    range = gt_range_offset(&range, 1); /* 1-based */
  if (!gt_str_length(sequenceid) ||
      (gt_cstr_table_get(srf->used_seqids, gt_str_get(sequenceid)) &&
       offset == GT_UNDEF_ULONG)) {
    /* sequenceid is empty or exists already (and no offset has been parsed)
       -> make one up */
    GtStr *seqid;
    char *base;
    base = gt_basename(gth_input_get_genomic_filename(input, filenum));
    seqid = gt_str_new_cstr(base);
    gt_free(base);
    gt_str_append_char(seqid, '|');
    gt_str_append_ulong(seqid, seqnum + 1); /* 1-based */
    seqid_store_add(srf->seqid_store, filenum, seqnum, seqid, offset);
    gt_assert(!gt_cstr_table_get(srf->used_seqids, gt_str_get(seqid)));
    gt_cstr_table_add(srf->used_seqids, gt_str_get(seqid));
    sr = gt_region_node_new(seqid_store_get(srf->seqid_store, filenum, seqnum),
                            range.start, range.end);
    gt_hashmap_add(sequence_regions,
                   (void*) gt_cstr_table_get(srf->used_seqids,
                                             gt_str_get(seqid)),
                   sr);
    gt_str_delete(seqid);
  }
  else {
    /* sequenceid does not exists already (or an offset has been parsed)
       -> use this one */
    if (!gt_cstr_table_get(srf->used_seqids, gt_str_get(sequenceid))) {
      /* no sequence region with this id exists -> create one */
      gt_cstr_table_add(srf->used_seqids, gt_str_get(sequenceid));
      seqid_store_add(srf->seqid_store, filenum, seqnum, sequenceid, offset);
      sr = gt_region_node_new(seqid_store_get(srf->seqid_store, filenum,
                                              seqnum), range.start, range.end);
      gt_hashmap_add(sequence_regions,
                     (void*) gt_cstr_table_get(srf->used_seqids,
                                              gt_str_get(sequenceid)),
                     sr);
    }
    else {
      GtRange prev_range, new_range;
      /* sequence region with this id exists already -> modify range */
      sr = gt_hashmap_get(sequence_regions, gt_str_get(sequenceid));
      gt_assert(sr);
      prev_range = gt_genome_node_get_range(sr);
      new_range = gt_range_join(&prev_range, &range);
      gt_genome_node_set_range(sr, &new_range);
      seqid_store_add(srf->seqid_store, filenum, seqnum, sequenceid, offset);
    }
  }
  gt_assert(sr);
}

static int accept_visitor(GT_UNUSED void *key, void *value, void *data,
                          GT_UNUSED GtError *err)
{
  GtGenomeNode *sr = value;
  GtNodeVisitor *visitor = data;
  int had_err;
  gt_error_check(err);
  gt_assert(sr && visitor);
  had_err = gt_genome_node_accept(sr, visitor, NULL);
  gt_assert(!had_err); /* should not happen */
  return 0;
}

void gth_region_factory_make(GthRegionFactory *srf, GtNodeVisitor *visitor,
                             GthInput *input)
{
  GtHashmap *sequence_regions;
  unsigned long i, j;
  GtStr *sequenceid;
  int had_err;
  gt_assert(srf && visitor && input);
  gt_assert(!srf->factory_was_used);
  srf->seqid_store = seqid_store_new(input);
  sequence_regions = gt_hashmap_new(GT_HASH_STRING, NULL,
                                    (GtFree) gt_genome_node_delete);
  sequenceid = gt_str_new();
  for (i = 0; i < gth_input_num_of_gen_files(input); i++) {
    for (j = 0; j < gth_input_num_of_gen_seqs(input, i); j++) {
      gt_str_reset(sequenceid);
      gth_input_save_gen_id(input, sequenceid, i, j);
      make_sequence_region(sequence_regions, sequenceid, srf, input, i, j);
    }
  }
  gt_str_delete(sequenceid);
  had_err = gt_hashmap_foreach_in_key_order(sequence_regions, accept_visitor,
                                            visitor, NULL);
  gt_assert(!had_err); /* should not happen */
  gt_hashmap_delete(sequence_regions);
  srf->factory_was_used = true;
}

GtStr* gth_region_factory_get_seqid(GthRegionFactory *srf,
                                    unsigned long filenum, unsigned long seqnum)
{
  gt_assert(srf && srf->factory_was_used);
  return seqid_store_get(srf->seqid_store, filenum, seqnum);
}

long gth_region_factory_offset(GthRegionFactory *srf,
                               unsigned long filenum, unsigned long seqnum)
{
  gt_assert(srf && srf->factory_was_used);
  return seqid_store_offset(srf->seqid_store, filenum, seqnum);
}
