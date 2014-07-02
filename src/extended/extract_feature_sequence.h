/*
  Copyright (c) 2006-2011 Gordon Gremme <gordon@gremme.org>
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

#ifndef EXTRACT_FEATURE_SEQUENCE_H
#define EXTRACT_FEATURE_SEQUENCE_H

#include "core/trans_table_api.h"
#include "extended/feature_node_api.h"
#include "extended/genome_node_api.h"
#include "extended/region_mapping_api.h"

int gt_extract_feature_sequence(GtStr *sequence, GtGenomeNode*,
                                const char *type, bool join, GtStr *seqid,
                                GtStrArray *target_ids, GtRegionMapping*,
                                GtError*);

/* Translates the sequences of all sub-features of type <type> appearing as
   direct children of <feature_node>. If <ttable> is not NULL, it will be used
   as the translation table. If <translation_fr1>, <translation_fr2> or
   <translation_fr3> are not NULL, they will contain the translation in the
   corresponding reading frame. Returns 0 on success, -1 otherwise and sets
   <err> accordingly. */
int gt_extract_and_translate_feature_sequence(GtFeatureNode *feature_node,
                                              const char *type, bool join,
                                              GtRegionMapping *rm,
                                              GtTransTable *ttable,
                                              GtStr *translation_fr1,
                                              GtStr *translation_fr2,
                                              GtStr *translation_fr3,
                                              GtError *err);

#endif
