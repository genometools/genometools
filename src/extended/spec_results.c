/*
  Copyright (c) 2014 Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2014 Genome Research Ltd.

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
#include "core/log_api.h"
#include "core/ma.h"
#include "core/str.h"
#include "core/str_array_api.h"
#include "core/symbol_api.h"
#include "core/unused_api.h"
#include "extended/comment_node_api.h"
#include "extended/feature_node_api.h"
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
  GtStrArray *failure_messages;
  GtStrArray *runtime_error_messages;
} GtSpecAspectNodeResult;

typedef struct {
  GtUword *succ, *fail, *err;
} GtSpecAspectCountInfo;

typedef struct {
  GtHashmap *aspect_node_results;
  GtStr *name;
  GtUword node_successes,
          node_failures,
          node_runtime_errors,
          successes,
          failures,
          runtime_errors;
} GtSpecAspect;

typedef struct {
  GtFile *outfile;
  bool details,
       colored,
       show_per_node;
  GtUword node_i;
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
  GtStrArray *warnings;
};

static GtSpecAspectNodeResult* gt_spec_aspect_node_result_new()
{
  GtSpecAspectNodeResult *sanr;

  sanr = gt_calloc(1, sizeof (*sanr));
  sanr->failure_messages = gt_str_array_new();
  sanr->runtime_error_messages = gt_str_array_new();
  return sanr;
}

static void gt_spec_aspect_node_result_delete(GtSpecAspectNodeResult *sanr)
{
  if (!sanr) return;
  gt_str_array_delete(sanr->failure_messages);
  gt_str_array_delete(sanr->runtime_error_messages);
  gt_free(sanr);
}

static GtSpecAspect* gt_spec_aspect_new(const char *name)
{
  GtSpecAspect *sa;
  gt_assert(name);

  sa = gt_calloc(1, sizeof (*sa));
  sa->name = gt_str_new_cstr(name);
  sa->aspect_node_results= gt_hashmap_new(GT_HASH_DIRECT,
                                    (GtFree) gt_genome_node_delete,
                                    (GtFree) gt_spec_aspect_node_result_delete);

  return sa;
}

static int gt_spec_aspect_report_node(void *key, void *value, void *data,
                                      GT_UNUSED GtError *err)
{
  GtGenomeNode *gn = (GtGenomeNode*) key;
  GtSpecAspectNodeResult *sanr = (GtSpecAspectNodeResult*) value;
  GtSpecResultsReportInfo *info = (GtSpecResultsReportInfo*) data;
  char tmpbuf[BUFSIZ];
  bool has_id = false;
  GtUword i = 0;

  if (gt_feature_node_try_cast(gn)) {
    const char *tmp;
    if ((tmp = gt_feature_node_get_attribute((GtFeatureNode*) gn, "ID"))) {
      has_id = true;
      (void) snprintf(tmpbuf, BUFSIZ, "%s, ", tmp);
    }
  }
  gt_file_xprintf(info->outfile,
                  "%s      offending node #" GT_WU " "
                  "(%sfrom %s, line %u):%s\n",
                  info->colored ? GT_SPEC_ANSI_COLOR_RED : "",
                  info->node_i++,
                  (has_id ? tmpbuf : ""),
                  gt_genome_node_get_filename(gn),
                  gt_genome_node_get_line_number(gn),
                  info->colored ? GT_SPEC_ANSI_COLOR_RESET : "");

  for (i = 0; i < gt_str_array_size(sanr->failure_messages); i++) {
    gt_file_xprintf(info->outfile,
                    "%s         %s%s\n",
                    info->colored ? GT_SPEC_ANSI_COLOR_RED : "",
                    gt_str_array_get(sanr->failure_messages, i),
                    info->colored ? GT_SPEC_ANSI_COLOR_RESET : "");
  }
  for (i = 0; i < gt_str_array_size(sanr->runtime_error_messages); i++) {
    gt_file_xprintf(info->outfile,
                    "%s         %s%s\n",
                    info->colored ? GT_SPEC_ANSI_COLOR_MAGENTA : "",
                    gt_str_array_get(sanr->runtime_error_messages, i),
                    info->colored ? GT_SPEC_ANSI_COLOR_RESET : "");
  }

  return 0;
}

static int gt_spec_aspect_count_stats(GT_UNUSED void *key, void *value,
                                      void *data, GT_UNUSED GtError *err)
{
  GtSpecAspectNodeResult *sanr = (GtSpecAspectNodeResult*) value;
  GtSpecAspectCountInfo *info = (GtSpecAspectCountInfo*) data;

  if (gt_str_array_size(sanr->runtime_error_messages) > 0) {
    (*info->err)++;
  } else if (gt_str_array_size(sanr->failure_messages) > 0) {
    (*info->fail)++;
  } else {
    (*info->succ)++;
  }
  return 0;
}

static void gt_spec_aspect_report(GtSpecAspect *sa, GtFile *outfile,
                                  GtSpecResultsReportInfo *info)
{
  GtSpecAspectCountInfo cinfo;
  GT_UNUSED int rval;
  gt_assert(sa && info);

  /* count per-node results */
  cinfo.succ = &(sa->node_successes);
  cinfo.fail = &(sa->node_failures);
  cinfo.err = &(sa->node_runtime_errors);
  rval = gt_hashmap_foreach(sa->aspect_node_results, gt_spec_aspect_count_stats,
                            &cinfo, NULL);
  gt_assert(rval == 0);

  if (info->show_per_node) {
    sa->successes = sa->node_successes;
    sa->failures = sa->node_failures;
    sa->runtime_errors = sa->node_runtime_errors;
  }

  /* output stats */
  info->node_i = 1;
  GtStr *outbuf = gt_str_new_cstr("  - ");
  gt_str_append_str(outbuf, sa->name);
  gt_str_append_cstr(outbuf, " (");
  if (sa->successes > 0) {
    if (info->colored)
      gt_str_append_cstr(outbuf, GT_SPEC_ANSI_COLOR_GREEN);
    gt_str_append_ulong(outbuf, sa->successes);
    gt_str_append_cstr(outbuf, " success");
    if (sa->successes > 1)
      gt_str_append_cstr(outbuf, "es");
    if (info->colored)
      gt_str_append_cstr(outbuf, GT_SPEC_ANSI_COLOR_RESET);
  }
  if (sa->failures > 0) {
    if (sa->successes > 0)
      gt_str_append_cstr(outbuf, ", ");
    if (info->colored)
      gt_str_append_cstr(outbuf, GT_SPEC_ANSI_COLOR_RED);
    gt_str_append_ulong(outbuf, sa->failures);
    gt_str_append_cstr(outbuf, " failure");
    if (sa->failures > 1)
      gt_str_append_cstr(outbuf, "s");
    if (info->colored)
      gt_str_append_cstr(outbuf, GT_SPEC_ANSI_COLOR_RESET);
  }
  if (sa->runtime_errors > 0) {
    if (sa->successes > 0 || sa->failures > 0)
      gt_str_append_cstr(outbuf, ", ");
    if (info->colored)
      gt_str_append_cstr(outbuf, GT_SPEC_ANSI_COLOR_MAGENTA);
    gt_str_append_ulong(outbuf, sa->runtime_errors);
    gt_str_append_cstr(outbuf, " runtime error");
    if (sa->runtime_errors > 1)
      gt_str_append_cstr(outbuf, "s");
    if (info->colored)
      gt_str_append_cstr(outbuf, GT_SPEC_ANSI_COLOR_RESET);
  }
  gt_str_append_cstr(outbuf, ")\n");
  gt_file_xprintf(outfile, "%s", gt_str_get(outbuf));
  gt_str_delete(outbuf);
  if (info->details) {
    rval = gt_hashmap_foreach(sa->aspect_node_results,
                              gt_spec_aspect_report_node, info, NULL);
    gt_assert(rval == 0);
  }
}

static void gt_spec_aspect_delete(GtSpecAspect *sa)
{
  if (!sa) return;
  gt_str_delete(sa->name);
  gt_hashmap_delete(sa->aspect_node_results);
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
  spec_results->warnings = gt_str_array_new();
  return spec_results;
}

void gt_spec_results_add_cc(GtSpecResults *sr)
{
  gt_assert(sr);
  sr->checked_ccs++;
}

void gt_spec_results_record_warning(GtSpecResults *sr, const char *w)
{
  gt_assert(sr && w);
  gt_str_array_add_cstr(sr->warnings, w);
}

void gt_spec_results_add_result(GtSpecResults *sr,
                                const char *aspect,
                                GtGenomeNode *node,
                                GtSpecResultStatus status,
                                const char *error_string)
{
  GtSpecAspect *sa = NULL;
  GtSpecAspectNodeResult *sanr = NULL;
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
  if (!(sanr = gt_hashmap_get(sa->aspect_node_results, node))) {
    sanr = gt_spec_aspect_node_result_new();
    node = gt_genome_node_ref(node);
    gt_hashmap_add(sa->aspect_node_results, node, sanr);
  }
  switch (status) {
    case GT_SPEC_SUCCESS:
      sa->successes++;
      break;
    case GT_SPEC_FAILURE:
      sa->failures++;
      gt_assert(sanr);
      gt_str_array_add_cstr(sanr->failure_messages, error_string);
      break;
    case GT_SPEC_RUNTIME_ERROR:
      sa->runtime_errors++;
      gt_assert(sanr);
      gt_str_array_add_cstr(sanr->runtime_error_messages, error_string);
      break;
  }
  sr->checked_nodes++;
}

static int gt_spec_results_report_single(GT_UNUSED void *key, void *value,
                                         void *data, GT_UNUSED GtError *err)
{
  GtSpecResultsReportInfo *info = (GtSpecResultsReportInfo*) data;
  gt_spec_aspect_report((GtSpecAspect*) value, info->outfile, info);
  return 0;
}

static int gt_spec_results_report_features(void *key, void *value, void *data,
                                           GT_UNUSED GtError *err)
{
  const char *type = (const char*) key;
  GtSpecResultsReportInfo *info = (GtSpecResultsReportInfo*) data;
  gt_file_xprintf(info->outfile, "a %s%s%s feature\n",
                           info->colored ? GT_SPEC_ANSI_COLOR_YELLOW : "",
                           type,
                           info->colored ? GT_SPEC_ANSI_COLOR_RESET : "");
  gt_hashmap_foreach((GtHashmap*) value, gt_spec_results_report_single,
                     info, NULL);
  return 0;
}

void gt_spec_results_report(GtSpecResults *sr, GtFile *outfile,
                            const char *specfile, bool details,bool colored,
                            bool show_per_node)
{
  GtSpecResultsReportInfo info;
  gt_assert(sr);

  info.outfile = outfile;
  info.details = details;
  info.colored = colored;
  info.show_per_node = show_per_node;

  if (sr->seen_feature || sr->seen_meta || sr->seen_region || sr->seen_comment
        || sr->seen_sequence) {
    gt_file_xprintf(outfile, "According to the specification in %s%s%s,\n\n",
                               colored ? GT_SPEC_ANSI_COLOR_YELLOW : "",
                               specfile,
                               colored ? GT_SPEC_ANSI_COLOR_RESET : "");
  }
  if (sr->seen_feature) {
    gt_hashmap_foreach_in_key_order(sr->feature_aspects,
                                    gt_spec_results_report_features,
                                    &info, NULL);
  }
  if (sr->seen_meta) {
    gt_file_xprintf(outfile, "a %smeta%s node\n",
                             colored ? GT_SPEC_ANSI_COLOR_YELLOW : "",
                             colored ? GT_SPEC_ANSI_COLOR_RESET : "");
    gt_hashmap_foreach_in_key_order(sr->meta_aspects,
                                    gt_spec_results_report_single,
                                    &info, NULL);
  }
  if (sr->seen_region) {
    gt_file_xprintf(outfile, "a %sregion%s node\n",
                             colored ? GT_SPEC_ANSI_COLOR_YELLOW : "",
                             colored ? GT_SPEC_ANSI_COLOR_RESET : "");
    gt_hashmap_foreach_in_key_order(sr->region_aspects,
                                    gt_spec_results_report_single,
                                    &info, NULL);
  }
  if (sr->seen_comment) {
    gt_file_xprintf(outfile, "a %scomment%s node\n",
                             colored ? GT_SPEC_ANSI_COLOR_YELLOW : "",
                             colored ? GT_SPEC_ANSI_COLOR_RESET : "");
    gt_hashmap_foreach_in_key_order(sr->comment_aspects,
                                    gt_spec_results_report_single,
                                    &info, NULL);
  }
  if (sr->seen_sequence) {
    gt_file_xprintf(outfile, "a %ssequence%s node\n",
                             colored ? GT_SPEC_ANSI_COLOR_YELLOW : "",
                             colored ? GT_SPEC_ANSI_COLOR_RESET : "");
    gt_hashmap_foreach_in_key_order(sr->sequence_aspects,
                                    gt_spec_results_report_single,
                                    &info, NULL);
  }
  if (gt_str_array_size(sr->warnings) > 0) {
    gt_file_xprintf(outfile, "\nThere ha%s been %s"GT_WU"%s parser warning%s\n",
                             gt_str_array_size(sr->warnings) > 1 ? "ve" : "s",
                             colored ? GT_SPEC_ANSI_COLOR_YELLOW : "",
                             gt_str_array_size(sr->warnings),
                             colored ? GT_SPEC_ANSI_COLOR_RESET : "",
                             gt_str_array_size(sr->warnings) > 1 ? "s." : ".");
    if (details) {
      GtUword i;
      for (i = 0; i < gt_str_array_size(sr->warnings); i++) {
        gt_file_xprintf(outfile, "  - %s\n", gt_str_array_get(sr->warnings, i));
      }

    }
  }
  gt_file_xprintf(outfile, "\nTraversed " GT_WU " CCs "
                           "(" GT_WU " feature types), "
                           "checked " GT_WU " nodes for " GT_WU " aspects.\n",
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
  gt_str_array_delete(spec_results->warnings);
  gt_free(spec_results);
}
