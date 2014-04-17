/*
  Copyright (c) 2014 Sascha Steinbiss <sascha@steinbiss.name>

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

#include <stdio.h>
#include "core/array.h"
#include "core/cstr_api.h"
#include "core/file.h"
#include "core/hashmap_api.h"
#include "core/ma.h"
#include "core/str.h"
#include "core/symbol_api.h"
#include "core/unused_api.h"
#include "extended/comment_node_api.h"
#include "extended/feature_node_api.h"
#include "extended/genome_node.h"
#include "extended/genome_node.h"
#include "extended/genome_node.h"
#include "extended/spec_results.h"

#define GT_SPEC_ANSI_COLOR_RED     "\x1b[31m"
#define GT_SPEC_ANSI_COLOR_GREEN   "\x1b[32m"
#define GT_SPEC_ANSI_COLOR_YELLOW  "\x1b[33m"
#define GT_SPEC_ANSI_COLOR_BLUE    "\x1b[34m"
#define GT_SPEC_ANSI_COLOR_MAGENTA "\x1b[35m"
#define GT_SPEC_ANSI_COLOR_CYAN    "\x1b[36m"
#define GT_SPEC_ANSI_COLOR_RESET   "\x1b[0m"

typedef struct {
  GtArray *failed_nodes;
  GtArray *failure_messages;
  GtStr *name;
  GtUword successes,
          failures;
} GtSpecAspect;

typedef struct {
  GtFile *outfile;
  bool details,
       colored;
} GtSpecResultsReportInfo;

struct GtSpecResults
{
  GtHashmap *feature_aspects,
            *meta_aspects,
            *region_aspects,
            *comment_aspects,
            *sequence_aspects;
  bool seen_feature,
       seen_meta,
       seen_region,
       seen_comment,
       seen_sequence;
  GtUword checked_types,
          checked_aspects,
          checked_ccs,
          checked_nodes;
};

static GtSpecAspect* gt_spec_aspect_new(const char *name)
{
  GtSpecAspect *sa;
  gt_assert(name);

  sa = gt_calloc(1, sizeof (*sa));
  sa->name = gt_str_new_cstr(name);
  sa->failed_nodes = gt_array_new(sizeof (GtGenomeNode*));
  sa->failure_messages = gt_array_new(sizeof (GtStr*));
  sa->successes = sa->failures = 0;

  return sa;
}

static void gt_spec_aspect_report(GtSpecAspect *sa, GtFile *outfile,
                                  bool details, bool colored)
{
  GtUword i;
  gt_assert(sa);
  if (gt_array_size(sa->failed_nodes) > 0) {
    gt_file_xprintf(outfile, "%s  - %s (" GT_WU " failed)%s\n",
                             colored ? GT_SPEC_ANSI_COLOR_RED : "",
                             gt_str_get(sa->name),
                             gt_array_size(sa->failed_nodes),
                             colored ? GT_SPEC_ANSI_COLOR_RESET : "");
    if (details) {
      for (i = 0; i < gt_array_size(sa->failed_nodes); i++) {
        char tmpbuf[BUFSIZ];
        bool has_id = false;
        GtGenomeNode *gn = *(GtGenomeNode**) gt_array_get(sa->failed_nodes, i);
        if (gt_feature_node_try_cast(gn)) {
          const char *tmp;
          if ((tmp = gt_feature_node_get_attribute((GtFeatureNode*) gn,
                                                   "ID"))) {
            has_id = true;
          }
          if (has_id)
            (void) snprintf(tmpbuf, BUFSIZ, "%s, ", tmp);
        }
        GtStr *errmsg = *(GtStr**) gt_array_get(sa->failure_messages, i);
        gt_file_xprintf(outfile,
                       "%s"
                       "      offending node #" GT_WU " (%sfrom %s, line %u):\n"
                       "         %s%s\n",
                       colored ? GT_SPEC_ANSI_COLOR_RED : "",
                       i + 1,
                       (has_id ? tmpbuf : ""),
                       gt_genome_node_get_filename(gn),
                       gt_genome_node_get_line_number(gn),
                       gt_str_get(errmsg),
                       colored ? GT_SPEC_ANSI_COLOR_RESET : "");
      }
    }
  } else {
    gt_file_xprintf(outfile, "%s  - %s%s\n",
                    colored ? GT_SPEC_ANSI_COLOR_GREEN : "",
                    gt_str_get(sa->name),
                    colored ? GT_SPEC_ANSI_COLOR_RESET : "");
  }
}

static void gt_spec_aspect_delete(GtSpecAspect *sa)
{
  GtUword i;
  if (!sa) return;
  gt_str_delete(sa->name);
  for (i = 0; i < gt_array_size(sa->failed_nodes); i++) {
    gt_genome_node_delete(*(GtGenomeNode**) gt_array_get(sa->failed_nodes, i));
  }
  gt_array_delete(sa->failed_nodes);
  for (i = 0; i < gt_array_size(sa->failure_messages); i++) {
    gt_str_delete(*(GtStr**) gt_array_get(sa->failure_messages, i));
  }
  gt_array_delete(sa->failure_messages);
  gt_free(sa);
}

GtSpecResults *gt_spec_results_new(void)
{
  GtSpecResults *spec_results;
  spec_results = gt_calloc(1, sizeof (GtSpecResults));
  spec_results->feature_aspects = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                                 (GtFree) gt_hashmap_delete);
  spec_results->meta_aspects = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                              (GtFree) gt_spec_aspect_delete);
  spec_results->sequence_aspects = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                               (GtFree) gt_spec_aspect_delete);
  spec_results->region_aspects = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                               (GtFree) gt_spec_aspect_delete);
  spec_results->comment_aspects = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                               (GtFree) gt_spec_aspect_delete);
  return spec_results;
}

void gt_spec_results_add_cc(GtSpecResults *sr)
{
  gt_assert(sr);
  sr->checked_ccs++;
}

void gt_spec_results_add_result(GtSpecResults *sr,
                                const char *aspect,
                                GtGenomeNode *node,
                                bool success,
                                const char *error_string)
{
  GtSpecAspect *sa = NULL;
  GtStr *errmsg;
  gt_assert(aspect && node && error_string);
  if (gt_feature_node_try_cast(node)) {
    GtHashmap *per_type_aspects;
    const char *type = gt_feature_node_get_type((GtFeatureNode*) node);
    if (!(per_type_aspects = gt_hashmap_get(sr->feature_aspects, type))) {
      per_type_aspects = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                        (GtFree) gt_spec_aspect_delete);
      gt_hashmap_add(sr->feature_aspects, gt_cstr_dup(type), per_type_aspects);
      sr->checked_types++;
      sr->seen_feature = true;
    }
    gt_assert(per_type_aspects);
    if (!(sa = gt_hashmap_get(per_type_aspects, aspect))) {
      sa = gt_spec_aspect_new(aspect);
      gt_hashmap_add(per_type_aspects, gt_cstr_dup(aspect), sa);
      sr->checked_aspects++;
    }
  } else if (gt_meta_node_try_cast(node)) {
    if (!(sa = gt_hashmap_get(sr->meta_aspects, aspect))) {
      sa = gt_spec_aspect_new(aspect);
      gt_hashmap_add(sr->meta_aspects, gt_cstr_dup(aspect), sa);
      sr->checked_aspects++;
      sr->seen_meta = true;
    }
  } else if (gt_region_node_try_cast(node)) {
    if (!(sa = gt_hashmap_get(sr->region_aspects, aspect))) {
      sa = gt_spec_aspect_new(aspect);
      gt_hashmap_add(sr->region_aspects, gt_cstr_dup(aspect), sa);
      sr->checked_aspects++;
      sr->seen_region = true;
    }
  } else if (gt_comment_node_try_cast(node)) {
    if (!(sa = gt_hashmap_get(sr->comment_aspects, aspect))) {
      sa = gt_spec_aspect_new(aspect);
      gt_hashmap_add(sr->comment_aspects, gt_cstr_dup(aspect), sa);
      sr->checked_aspects++;
      sr->seen_comment = true;
    }
  } else if (gt_sequence_node_try_cast(node)) {
    if (!(sa = gt_hashmap_get(sr->sequence_aspects, aspect))) {
      sa = gt_spec_aspect_new(aspect);
      gt_hashmap_add(sr->sequence_aspects, gt_cstr_dup(aspect), sa);
      sr->checked_aspects++;
      sr->seen_sequence = true;
    }
  }
  gt_assert(sa);
  if (success)
    sa->successes++;
  else {
    sa->failures++;
    node = gt_genome_node_ref(node);
    gt_array_add(sa->failed_nodes, node);
    errmsg = gt_str_new_cstr(error_string);
    gt_array_add(sa->failure_messages, errmsg);
  }
  sr->checked_nodes++;
}

static int gt_spec_results_report_single(GT_UNUSED void *key, void *value,
                                         void *data, GT_UNUSED GtError *err)
{
  GtSpecResultsReportInfo *info = (GtSpecResultsReportInfo*) data;
  gt_spec_aspect_report((GtSpecAspect*) value, info->outfile, info->details,
                        info->colored);
  return 0;
}

static int gt_spec_results_report_features(void *key, void *value, void *data,
                                           GT_UNUSED GtError *err)
{
  const char *type = (const char*) key;
  GtSpecResultsReportInfo *info = (GtSpecResultsReportInfo*) data;
  gt_file_xprintf(info->outfile, "Features of type %s%s%s\n",
                           info->colored ? GT_SPEC_ANSI_COLOR_YELLOW : "",
                           type,
                           info->colored ? GT_SPEC_ANSI_COLOR_RESET : "");
  gt_hashmap_foreach((GtHashmap*) value, gt_spec_results_report_single,
                     info, NULL);
  return 0;
}

void gt_spec_results_report(GtSpecResults *sr, GtFile *outfile, bool details,
                            bool colored)
{
  GtSpecResultsReportInfo info;
  gt_assert(sr);
  info.outfile = outfile;
  info.details = details;
  info.colored = colored;
  if (sr->seen_feature) {
    gt_hashmap_foreach_in_key_order(sr->feature_aspects,
                                    gt_spec_results_report_features,
                                    &info, NULL);
  }
  if (sr->seen_meta) {
    gt_file_xprintf(outfile, "%sMeta%s nodes\n",
                             colored ? GT_SPEC_ANSI_COLOR_YELLOW : "",
                             colored ? GT_SPEC_ANSI_COLOR_RESET : "");
    gt_hashmap_foreach_in_key_order(sr->meta_aspects,
                                    gt_spec_results_report_single,
                                    &info, NULL);
  }
  if (sr->seen_region) {
    gt_file_xprintf(outfile, "%sRegion%s nodes\n",
                             colored ? GT_SPEC_ANSI_COLOR_YELLOW : "",
                             colored ? GT_SPEC_ANSI_COLOR_RESET : "");
    gt_hashmap_foreach_in_key_order(sr->region_aspects,
                                    gt_spec_results_report_single,
                                    &info, NULL);
  }
  if (sr->seen_comment) {
    gt_file_xprintf(outfile, "%sComment%s nodes\n",
                             colored ? GT_SPEC_ANSI_COLOR_YELLOW : "",
                             colored ? GT_SPEC_ANSI_COLOR_RESET : "");
    gt_hashmap_foreach_in_key_order(sr->comment_aspects,
                                    gt_spec_results_report_single,
                                    &info, NULL);
  }
  if (sr->seen_sequence) {
    gt_file_xprintf(outfile, "%sSequence%s nodes\n",
                             colored ? GT_SPEC_ANSI_COLOR_YELLOW : "",
                             colored ? GT_SPEC_ANSI_COLOR_RESET : "");
    gt_hashmap_foreach_in_key_order(sr->sequence_aspects,
                                    gt_spec_results_report_single,
                                    &info, NULL);
  }
  gt_file_xprintf(outfile, "\nChecked " GT_WU " CCs (" GT_WU " feature types), "
                           GT_WU " nodes, " GT_WU " aspects.\n",
                           sr->checked_ccs, sr->checked_types,
                           sr->checked_nodes, sr->checked_aspects);
}

void gt_spec_results_delete(GtSpecResults *spec_results)
{
  if (!spec_results) return;
  gt_hashmap_delete(spec_results->feature_aspects);
  gt_hashmap_delete(spec_results->meta_aspects);
  gt_hashmap_delete(spec_results->comment_aspects);
  gt_hashmap_delete(spec_results->sequence_aspects);
  gt_hashmap_delete(spec_results->region_aspects);
  gt_free(spec_results);
}
