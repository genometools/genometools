/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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
#include "libgtext/regionmapping.h"

struct RegionMapping {
  Str *sequence_filename,
      *sequence_file, /* the (current) sequence file */
      *sequence_name; /* the (current) sequence name */
  Mapping *mapping;
  Bioseq *bioseq; /* the current bioseq */
};

RegionMapping* regionmapping_new_mapping(Str *mapping_filename, Env *env)
{
  RegionMapping *rm;
  env_error_check(env);
  assert(mapping_filename);
  rm = ma_calloc(1, sizeof (RegionMapping));
  rm->mapping = mapping_new(mapping_filename, "mapping", MAPPINGTYPE_STRING,
                            env);
  if (!rm->mapping) {
    regionmapping_delete(rm, env);
    return NULL;
  }
  return rm;
}

RegionMapping* regionmapping_new_seqfile(Str *sequence_filename, Env *env)
{
  RegionMapping *rm;
  assert(sequence_filename && env);
  rm = ma_calloc(1, sizeof (RegionMapping));
  rm->sequence_filename = str_ref(sequence_filename);
  return rm;
}

static Str* regionmapping_map(RegionMapping *rm, const char *sequence_region,
                              Env *env)
{
  env_error_check(env);
  assert(rm && sequence_region);
  if (rm->sequence_filename)
    return str_ref(rm->sequence_filename);
  else
    return mapping_map_string(rm->mapping, sequence_region, env);
}

static int update_bioseq_if_necessary(RegionMapping *rm, Str *seqid, Env *env)
{
  int had_err = 0;
  env_error_check(env);
  assert(rm && seqid);
  if (!rm->sequence_file || str_cmp(rm->sequence_name, seqid)) {
    str_delete(rm->sequence_file);
    rm->sequence_file = regionmapping_map(rm, str_get(seqid), env);
    if (!rm->sequence_file)
      had_err = -1;
    else {
      if (!rm->sequence_name)
        rm->sequence_name = str_new();
      else
        str_reset(rm->sequence_name);
      str_append_str(rm->sequence_name, seqid);
      bioseq_delete(rm->bioseq, env);
      rm->bioseq = bioseq_new_str(rm->sequence_file, env);
      if (!rm->bioseq)
        had_err = -1;
    }
  }
  return had_err;
}

int regionmapping_get_raw_sequence(RegionMapping *rm, const char **raw,
                                   Str *seqid, Env *env)
{
  int had_err = 0;
  env_error_check(env);
  assert(rm && seqid);
  had_err = update_bioseq_if_necessary(rm, seqid, env);
  if (!had_err)
    *raw = bioseq_get_raw_sequence(rm->bioseq);
  return had_err;
}

int regionmapping_get_raw_sequence_length(RegionMapping *rm,
                                          unsigned long *length, Str *seqid,
                                          Env *env)
{
  int had_err = 0;
  env_error_check(env);
  assert(rm && seqid);
  had_err = update_bioseq_if_necessary(rm, seqid, env);
  if (!had_err)
    *length = bioseq_get_raw_sequence_length(rm->bioseq);
  return had_err;
}

void regionmapping_delete(RegionMapping *rm, Env *env)
{
  if (!rm) return;
  str_delete(rm->sequence_filename);
  str_delete(rm->sequence_file);
  str_delete(rm->sequence_name);
  mapping_delete(rm->mapping, env);
  bioseq_delete(rm->bioseq, env);
  ma_free(rm);
}
