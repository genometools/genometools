/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <assert.h>
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#include "libgtcore/bioseq.h"
#include "libgtcore/ma.h"
#include "libgtext/mapping.h"
#include "libgtext/region_mapping.h"

struct RegionMapping {
  Str *sequence_filename,
      *sequence_file, /* the (current) sequence file */
      *sequence_name; /* the (current) sequence name */
  Mapping *mapping;
  Bioseq *bioseq; /* the current bioseq */
  unsigned int reference_count;
};

RegionMapping* region_mapping_new_mapping(Str *mapping_filename, Error *e)
{
  RegionMapping *rm;
  error_check(e);
  assert(mapping_filename);
  rm = ma_calloc(1, sizeof (RegionMapping));
  rm->mapping = mapping_new(mapping_filename, "mapping", MAPPINGTYPE_STRING, e);
  if (!rm->mapping) {
    region_mapping_delete(rm);
    return NULL;
  }
  return rm;
}

RegionMapping* region_mapping_new_seqfile(Str *sequence_filename)
{
  RegionMapping *rm;
  assert(sequence_filename);
  rm = ma_calloc(1, sizeof (RegionMapping));
  rm->sequence_filename = str_ref(sequence_filename);
  return rm;
}

RegionMapping* region_mapping_ref(RegionMapping *rm)
{
  assert(rm);
  rm->reference_count++;
  return rm;
}

static Str* region_mapping_map(RegionMapping *rm, const char *sequence_region,
                               Error *e)
{
  error_check(e);
  assert(rm && sequence_region);
  if (rm->sequence_filename)
    return str_ref(rm->sequence_filename);
  else
    return mapping_map_string(rm->mapping, sequence_region, e);
}

static int update_bioseq_if_necessary(RegionMapping *rm, Str *seqid, Error *e)
{
  int had_err = 0;
  error_check(e);
  assert(rm && seqid);
  if (!rm->sequence_file || str_cmp(rm->sequence_name, seqid)) {
    str_delete(rm->sequence_file);
    rm->sequence_file = region_mapping_map(rm, str_get(seqid), e);
    if (!rm->sequence_file)
      had_err = -1;
    else {
      if (!rm->sequence_name)
        rm->sequence_name = str_new();
      else
        str_reset(rm->sequence_name);
      str_append_str(rm->sequence_name, seqid);
      bioseq_delete(rm->bioseq);
      rm->bioseq = bioseq_new_str(rm->sequence_file, e);
      if (!rm->bioseq)
        had_err = -1;
    }
  }
  return had_err;
}

int region_mapping_get_raw_sequence(RegionMapping *rm, const char **raw,
                                    Str *seqid, Error *e)
{
  int had_err = 0;
  error_check(e);
  assert(rm && seqid);
  had_err = update_bioseq_if_necessary(rm, seqid, e);
  if (!had_err)
    *raw = bioseq_get_raw_sequence(rm->bioseq);
  return had_err;
}

int region_mapping_get_raw_sequence_length(RegionMapping *rm,
                                           unsigned long *length, Str *seqid,
                                           Error *e)
{
  int had_err = 0;
  error_check(e);
  assert(rm && seqid);
  had_err = update_bioseq_if_necessary(rm, seqid, e);
  if (!had_err)
    *length = bioseq_get_raw_sequence_length(rm->bioseq);
  return had_err;
}

void region_mapping_delete(RegionMapping *rm)
{
  if (!rm) return;
  if (rm->reference_count) {
    rm->reference_count--;
    return;
  }
  str_delete(rm->sequence_filename);
  str_delete(rm->sequence_file);
  str_delete(rm->sequence_name);
  mapping_delete(rm->mapping);
  bioseq_delete(rm->bioseq);
  ma_free(rm);
}
