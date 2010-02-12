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
#include "core/parseutils.h"
#include "core/undef.h"
#include "gth/sequence_region_factory.h"

typedef struct {
  unsigned long num_of_files,
                *num_of_sequences;
  GtStr ***store;
  long **offsets;
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
      ss->offsets[i][j] = GT_UNDEF_LONG;
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
                            unsigned long seqnum, GtStr *seqid, long offset)
{
  gt_assert(ss && seqid);
  gt_assert(gt_str_length(seqid)); /* is not empty */
  gt_assert(filenum < ss->num_of_files);
  gt_assert(seqnum < ss->num_of_sequences[filenum]);
  gt_assert(!ss->store[filenum][seqnum]); /* is unused */
  ss->store[filenum][seqnum] = gt_str_clone(seqid);
  ss->offsets[filenum][seqnum] = offset == GT_UNDEF_LONG ? 1 : offset;
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

static long seqid_store_offset(SeqidStore *ss, unsigned long filenum,
                               unsigned long seqnum)
{
  long offset;
  gt_assert(ss);
  gt_assert(filenum < ss->num_of_files);
  gt_assert(seqnum < ss->num_of_sequences[filenum]);
  offset = ss->offsets[filenum][seqnum];
  gt_assert(offset != GT_UNDEF_LONG); /* is defined */
  return offset;
}

struct SequenceRegionFactory{
  GtCstrTable *used_seqids;
  bool factory_was_used,
       use_desc_ranges;
  SeqidStore *seqid_store;
};

SequenceRegionFactory* sequence_region_factory_new(bool use_desc_ranges)
{
  SequenceRegionFactory *srf = gt_calloc(1, sizeof *srf);
  srf->used_seqids = gt_cstr_table_new();
  srf->use_desc_ranges = use_desc_ranges;
  return srf;
}

void sequence_region_factory_delete(SequenceRegionFactory *srf)
{
  if (!srf) return;
  gt_cstr_table_delete(srf->used_seqids);
  seqid_store_delete(srf->seqid_store);
  gt_free(srf);
}

/* Range descriptions have the folowing format: III:1000001..2000000
   That is, the part between ':' and '..' denotes the offset. */
static long parse_desc_range(GthInput *input,
                             unsigned long filenum, unsigned long seqnum)
{
  long offset;
  GtStr *description;
  unsigned long i;
  char *desc;
  gt_assert(input);
  description = gt_str_new();
  gth_input_get_genomic_description(input, description, filenum, seqnum);
  desc = gt_str_get(description);
  /* find ':' */
  for (i = 0; i < gt_str_length(description); i++) {
    if (desc[i] == ':')
      break;
  }
  if (i == gt_str_length(description)) {
    /* no ':' found */
    gt_str_delete(description);
    return GT_UNDEF_LONG;
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
    gt_str_delete(description);
    return GT_UNDEF_LONG;
  }
  /* parse range */
  gt_assert(desc[i-1] == '.' && desc[i] == '.');
  desc[i-1] = '\0';
  if (gt_parse_long(&offset, desc)) {
    /* parsing failed */
    gt_str_delete(description);
    return GT_UNDEF_LONG;
  }
  gt_str_delete(description);
  return offset;
}

static GtGenomeNode *make_sequence_region(GtStr *sequenceid,
                                          SequenceRegionFactory *srf,
                                          GthInput *input,
                                          unsigned long filenum,
                                          unsigned long seqnum)
{
  GtGenomeNode *sr = NULL;
  GtRange range;
  long offset = GT_UNDEF_LONG;
  gt_assert(sequenceid && srf && input);
  if (gth_input_use_substring_spec(input)) {
    range.start = gth_input_genomic_substring_from(input);
    range.end   = gth_input_genomic_substring_to(input);
  }
  else {
    range = gth_input_get_relative_genomic_range(input, filenum, seqnum);
  }
  if (srf->use_desc_ranges)
    offset = parse_desc_range(input, filenum, seqnum);
  if (offset != GT_UNDEF_LONG)
    range = gt_range_offset(&range, offset);
  else
    range = gt_range_offset(&range, 1); /* 1-based */
  if (!gt_str_length(sequenceid) ||
      (gt_cstr_table_get(srf->used_seqids, gt_str_get(sequenceid)) &&
       offset == GT_UNDEF_LONG)) {
    /* sequenceid is empty or exists already (and no offset has been parsed)
       -> make one up */
    GtStr *seqid;
    char *base;
    base = gt_basename(gth_input_get_genomic_filename(input, filenum));
    seqid = gt_str_new_cstr(base);
    gt_free(base);
    gt_str_append_char(seqid, '|');
    gt_str_append_ulong(seqid, seqnum + 1); /* 1-based */
    sr = gt_region_node_new(seqid, range.start, range.end);
    seqid_store_add(srf->seqid_store, filenum, seqnum, seqid, offset);
    gt_str_delete(seqid);
  }
  else {
    /* sequenceid does not exists already (or an offset has been parsed)
       -> use this one */
    sr = gt_region_node_new(sequenceid, range.start, range.end);
    seqid_store_add(srf->seqid_store, filenum, seqnum, sequenceid, offset);
    if (!gt_cstr_table_get(srf->used_seqids, gt_str_get(sequenceid)))
      gt_cstr_table_add(srf->used_seqids, gt_str_get(sequenceid));
  }
  gt_assert(sr);
  return sr;
}

void sequence_region_factory_make(SequenceRegionFactory *srf,
                                  GtNodeVisitor *visitor,
                                  GthInput *input)
{
  unsigned long i, j;
  GtStr *sequenceid;
  GtGenomeNode *sr;
  int had_err;
  gt_assert(srf && visitor && input);
  gt_assert(!srf->factory_was_used);
  srf->seqid_store = seqid_store_new(input);
  sequenceid = gt_str_new();
  for (i = 0; i < gth_input_num_of_gen_files(input); i++) {
    for (j = 0; j < gth_input_num_of_gen_seqs(input, i); j++) {
      gt_str_reset(sequenceid);
      gth_input_save_gen_id(input, sequenceid, i, j);
      sr = make_sequence_region(sequenceid, srf, input, i, j);
      had_err = gt_genome_node_accept(sr, visitor, NULL);
      gt_assert(!had_err); /* should not happen */
      gt_genome_node_delete(sr);
    }
  }
  gt_str_delete(sequenceid);
  srf->factory_was_used = true;
}

GtStr* sequence_region_factory_get_seqid(SequenceRegionFactory *srf,
                                         unsigned long filenum,
                                         unsigned long seqnum)
{
  gt_assert(srf && srf->factory_was_used);
  return seqid_store_get(srf->seqid_store, filenum, seqnum);
}

long sequence_region_factory_offset(SequenceRegionFactory *srf,
                                    unsigned long filenum,
                                    unsigned long seqnum)
{
  gt_assert(srf && srf->factory_was_used);
  return seqid_store_offset(srf->seqid_store, filenum, seqnum);
}
