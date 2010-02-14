/*
  Copyright (c) 2007-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#include "core/assert_api.h"
#include "core/bioseq.h"
#include "core/ma.h"
#include "extended/mapping.h"
#include "extended/region_mapping.h"

struct GtRegionMapping {
  GtStr *sequence_filename,
        *sequence_file, /* the (current) sequence file */
        *sequence_name; /* the (current) sequence name */
  bool usedesc;
  GtMapping *mapping;
  GtBioseq *bioseq; /* the current bioseq */
  unsigned int reference_count;
};

GtRegionMapping* gt_region_mapping_new_mapping(GtStr *mapping_filename,
                                               GtError *err)
{
  GtRegionMapping *rm;
  gt_error_check(err);
  gt_assert(mapping_filename);
  rm = gt_calloc(1, sizeof (GtRegionMapping));
  rm->mapping = gt_mapping_new(mapping_filename, "mapping", MAPPINGTYPE_STRING,
                               err);
  if (!rm->mapping) {
    gt_region_mapping_delete(rm);
    return NULL;
  }
  return rm;
}

GtRegionMapping* gt_region_mapping_new_seqfile(GtStr *sequence_filename,
                                               bool usedesc)
{
  GtRegionMapping *rm;
  gt_assert(sequence_filename);
  rm = gt_calloc(1, sizeof (GtRegionMapping));
  rm->sequence_filename = gt_str_ref(sequence_filename);
  rm->usedesc = usedesc;
  return rm;
}

GtRegionMapping* gt_region_mapping_ref(GtRegionMapping *rm)
{
  gt_assert(rm);
  rm->reference_count++;
  return rm;
}

static GtStr* region_mapping_map(GtRegionMapping *rm,
                                 const char *sequence_region, GtError *err)
{
  gt_error_check(err);
  gt_assert(rm && sequence_region);
  if (rm->sequence_filename)
    return gt_str_ref(rm->sequence_filename);
  else
    return gt_mapping_map_string(rm->mapping, sequence_region, err);
}

static int update_bioseq_if_necessary(GtRegionMapping *rm, GtStr *seqid,
                                      GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  gt_assert(rm && seqid);
  if (!rm->sequence_file || gt_str_cmp(rm->sequence_name, seqid)) {
    gt_str_delete(rm->sequence_file);
    rm->sequence_file = region_mapping_map(rm, gt_str_get(seqid), err);
    if (!rm->sequence_file)
      had_err = -1;
    else {
      if (!rm->sequence_name)
        rm->sequence_name = gt_str_new();
      else
        gt_str_reset(rm->sequence_name);
      gt_str_append_str(rm->sequence_name, seqid);
      gt_bioseq_delete(rm->bioseq);
      rm->bioseq = gt_bioseq_new_str(rm->sequence_file, err);
      if (!rm->bioseq)
        had_err = -1;
    }
  }
  return had_err;
}

int gt_region_mapping_get_raw_sequence(GtRegionMapping *rm, const char **rawseq,
                                       unsigned long *length, GtStr *seqid,
                                       GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  gt_assert(rm && rawseq && length && seqid);
  had_err = update_bioseq_if_necessary(rm, seqid, err);
  if (!had_err) {
    *rawseq = gt_bioseq_get_raw_sequence(rm->bioseq);
    *length = gt_bioseq_get_raw_sequence_length(rm->bioseq);
  }
  return had_err;
}

void gt_region_mapping_delete(GtRegionMapping *rm)
{
  if (!rm) return;
  if (rm->reference_count) {
    rm->reference_count--;
    return;
  }
  gt_str_delete(rm->sequence_filename);
  gt_str_delete(rm->sequence_file);
  gt_str_delete(rm->sequence_name);
  gt_mapping_delete(rm->mapping);
  gt_bioseq_delete(rm->bioseq);
  gt_free(rm);
}
