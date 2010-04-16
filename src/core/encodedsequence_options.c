/*
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#include "core/encodedsequence_options.h"
#include "core/ma.h"

struct GtEncodedsequenceOptions {
  GtProgressTimer *progress_timer;
  GtStr *indexname,
        *symbolmap,
        *accesstype;
  GtStrArray *infilenames;
  bool isdna,
       isprotein,
       isplain,
       tistable,
       sdstable,
       ssptable,
       destable,
       withrange;
  GtLogger *logger;
};

GtEncodedsequenceOptions* gt_encodedsequence_options_new(void)
{
  GtEncodedsequenceOptions *o = gt_calloc(1, sizeof (GtEncodedsequenceOptions));
  /* init everything to NULL or false */
  o->withrange = true;
  return o;
}

void gt_encodedsequence_options_set_progress_timer(GtEncodedsequenceOptions *o,
                                                   GtProgressTimer *pt)
{
  gt_assert(o);
  o->progress_timer = pt;
}

GtProgressTimer* gt_encodedsequence_options_get_progress_timer(
                                                   GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  return o->progress_timer;
}

void gt_encodedsequence_options_set_indexname(GtEncodedsequenceOptions *o,
                                              GtStr *indexname)
{
  gt_assert(o);
  if (o->indexname)
    gt_str_delete(o->indexname);
  o->indexname = indexname;
  if (indexname)
    o->indexname = gt_str_ref(o->indexname);
}

GtStr* gt_encodedsequence_options_get_indexname(GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  return o->indexname;
}

void gt_encodedsequence_options_set_symbolmap_file(GtEncodedsequenceOptions *o,
                                                   GtStr *smapfile)
{
  gt_assert(o);
  if (o->symbolmap)
    gt_str_delete(o->symbolmap);
  o->symbolmap = smapfile;
  if (smapfile)
    o->symbolmap = gt_str_ref(o->symbolmap);
}
GtStr* gt_encodedsequence_options_get_symbolmap_file(
                                                   GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  return o->symbolmap;
}

void gt_encodedsequence_options_set_access_type(GtEncodedsequenceOptions *o,
                                                GtStr *str_sat)
{
  gt_assert(o);
  if (o->accesstype)
    gt_str_delete(o->accesstype);
  o->accesstype = str_sat;
  if (str_sat)
    o->accesstype = gt_str_ref(o->accesstype);
}

GtStr* gt_encodedsequence_options_get_access_type(GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  return o->accesstype;
}

void gt_encodedsequence_options_set_input_sequences(GtEncodedsequenceOptions *o,
                                                    GtStrArray *filenametab)
{
  gt_assert(o);
  if (o->infilenames)
    gt_str_array_delete(o->infilenames);
  o->infilenames = filenametab;
  if (filenametab)
    o->infilenames = gt_str_array_ref(o->infilenames);
}
GtStrArray* gt_encodedsequence_options_get_input_sequences(
                                                   GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  return o->infilenames;
}

void gt_encodedsequence_options_set_input_dna(GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  o->isdna = true;
  o->isprotein = false;
  o->isplain = false;
}

bool gt_encodedsequence_options_get_input_dna(GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  return o->isdna;
}

void gt_encodedsequence_options_set_input_protein(GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  o->isdna = false;
  o->isprotein = true;
  o->isplain = false;
}

bool gt_encodedsequence_options_get_input_protein(GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  return o->isprotein;
}

void gt_encodedsequence_options_set_input_plain(GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  o->isdna = false;
  o->isprotein = false;
  o->isplain = true;
}

bool gt_encodedsequence_options_get_input_plain(GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  return o->isplain;
}

void gt_encodedsequence_options_enable_tis_table_usage(
                                                   GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  o->tistable = true;
}

void gt_encodedsequence_options_disable_tis_table_usage(
                                                   GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  o->tistable = false;
}

bool gt_encodedsequence_options_get_tis_table_usage(GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  return o->tistable;
}

void gt_encodedsequence_options_enable_des_table_usage(
                                                   GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  o->destable = true;
}

void gt_encodedsequence_options_disable_des_table_usage(
                                                   GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  o->destable = false;
}

bool gt_encodedsequence_options_get_des_table_usage(GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  return o->destable;
}

void gt_encodedsequence_options_enable_sds_table_usage(
                                                   GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  o->sdstable = true;
}

void gt_encodedsequence_options_disable_sds_table_usage(
                                                   GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  o->sdstable = false;
}

bool gt_encodedsequence_options_get_sds_table_usage(GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  return o->sdstable;
}

void gt_encodedsequence_options_enable_ssp_table_usage(
                                                   GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  o->ssptable = true;
}

void gt_encodedsequence_options_disable_ssp_table_usage(
                                                   GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  o->ssptable = false;
}

bool gt_encodedsequence_options_get_ssp_table_usage(GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  return o->ssptable;
}

void gt_encodedsequence_options_enable_range_iteration(
                                                   GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  o->withrange = true;
}

void gt_encodedsequence_options_disable_range_iteration(
                                                   GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  o->withrange = false;
}

bool gt_encodedsequence_options_get_range_iteration(GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  return o->withrange;
}

void gt_encodedsequence_options_set_logger(GtEncodedsequenceOptions *o,
                                           GtLogger *logger)
{
  gt_assert(o);
  o->logger = logger;
}

GtLogger* gt_encodedsequence_options_get_logger(GtEncodedsequenceOptions *o)
{
  gt_assert(o);
  return o->logger;
}

void gt_encodedsequence_options_delete(GtEncodedsequenceOptions *o)
{
  if (!o) return;
  if (o->indexname)
    gt_str_delete(o->indexname);
  if (o->symbolmap)
    gt_str_delete(o->symbolmap);
  if (o->accesstype)
    gt_str_delete(o->accesstype);
  if (o->infilenames)
    gt_str_array_delete(o->infilenames);
  gt_free(o);
}
