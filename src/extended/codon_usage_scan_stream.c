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

#include <ctype.h>
#include "core/array_api.h"
#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/queue.h"
#include "core/symbol_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/feature_node.h"
#include "extended/feature_type.h"
#include "extended/merge_stream.h"
#include "extended/codon_usage_scan_stream.h"

struct GtCodonUsageScanStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtArray *cur_cds_set;
  GtRange cur_cds_range;
  GtQueue *outqueue;
};

#define gt_codon_usage_scan_stream_cast(GS)\
        gt_node_stream_cast(gt_codon_usage_scan_stream_class(), GS);

static int codon_usage_scan_stream_process_current_cds(
                                               GtCodonUsageScanStream *cuss,
                                               GtError *err)
{
  int had_err = 0;
  GtUword i;
  double max = DBL_MIN;
  GtFeatureNode *maxnode = NULL;
  GtUword nof_cdss = gt_array_size(cuss->cur_cds_set);
  gt_error_check(err);

  if (nof_cdss > 0) {
    gt_assert(gt_queue_size(cuss->outqueue) == 0);
    for (i = 0; !had_err && i < nof_cdss; i++) {
      GtFeatureNode *cds;
      cds = *(GtFeatureNode**) gt_array_get(cuss->cur_cds_set, i);

      if (gt_feature_node_score_is_defined(cds)) {
        if (gt_feature_node_get_score(cds) > max) {
          max = gt_feature_node_get_score(cds);
          maxnode = cds;
        }
      }
    }
    if (maxnode) {
      gt_genome_node_ref((GtGenomeNode*) maxnode);
      gt_queue_add(cuss->outqueue, maxnode);
    }
  } else {
    /* no nodes for this cds cluster, delete it */
    for (i = 0; !had_err && i < nof_cdss; i++) {
      gt_genome_node_delete(*(GtGenomeNode**) gt_array_get(cuss->cur_cds_set,
                                                           i));
    }
  }
  gt_array_reset(cuss->cur_cds_set);
  return had_err;
}

static int codon_usage_scan_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                        GtError *err)
{
  GtCodonUsageScanStream *cuss;
  int had_err = 0;
  bool complete_cluster = false;
  GtGenomeNode *mygn = NULL;
  GtFeatureNode *fn = NULL;
  const char *cds = gt_symbol(gt_ft_CDS);
  gt_error_check(err);
  cuss = gt_codon_usage_scan_stream_cast(ns);

  /* if there are still nodes left in the buffer, output them */
  if (gt_queue_size(cuss->outqueue) > 0) {
    *gn = (GtGenomeNode*) gt_queue_get(cuss->outqueue);
    return had_err;
  } else complete_cluster = false;

  while (!had_err && !complete_cluster) {
    had_err = gt_node_stream_next(cuss->in_stream, &mygn, err);

    /* stop if stream is at the end */
    if (had_err || !mygn) break;

    /* process all feature nodes */
    if ((fn = gt_feature_node_try_cast(mygn))) {
      GtGenomeNode *addgn;
      const char *type = gt_feature_node_get_type(fn);
      GtRange new_rng = gt_genome_node_get_range(mygn);
      if (type == cds) {
        /* -----> this is a cds <----- */
        if (gt_array_size(cuss->cur_cds_set) == 0UL) {
          /* new overlapping cds cluster */
          addgn = gt_genome_node_ref(mygn);
          gt_array_add(cuss->cur_cds_set, addgn);
          cuss->cur_cds_range = gt_genome_node_get_range(mygn);
        } else {
          if (gt_range_overlap(&new_rng, &cuss->cur_cds_range)) {
            /* cds overlaps with current one, add to cluster */
            addgn = gt_genome_node_ref(mygn);
            gt_array_add(cuss->cur_cds_set, addgn);
            cuss->cur_cds_range = gt_range_join(&cuss->cur_cds_range,
                                                   &new_rng);
          } else {
            /* finish current cluster and start a new one */
            had_err = codon_usage_scan_stream_process_current_cds(cuss, err);
            if (!had_err) {
              addgn = gt_genome_node_ref(mygn);
              gt_array_add(cuss->cur_cds_set, addgn);
              cuss->cur_cds_range = gt_genome_node_get_range(mygn);
            }
            if (gt_queue_size(cuss->outqueue) > 0) {
              *gn = (GtGenomeNode*) gt_queue_get(cuss->outqueue);
              complete_cluster = true;
            }
          }
        }
        /* from now on, cdss are kept in cds cluster arrays only */
        gt_genome_node_delete(mygn);
      }
    } else {
      /* other nodes */
      had_err = codon_usage_scan_stream_process_current_cds(cuss, err);
      if (!had_err) {
        gt_queue_add(cuss->outqueue, mygn);
      }
      if (gt_queue_size(cuss->outqueue) > 0) {
        *gn = (GtGenomeNode*) gt_queue_get(cuss->outqueue);
        complete_cluster = true;
      }
    }
  }

  return had_err;
}

static void codon_usage_scan_stream_free(GtNodeStream *ns)
{
  GtUword i;
  GtCodonUsageScanStream *cuss;
  if (!ns) return;
  cuss = gt_codon_usage_scan_stream_cast(ns);
  while (gt_queue_size(cuss->outqueue) > 0) {
    gt_genome_node_delete((GtGenomeNode*) gt_queue_get(cuss->outqueue));
  }
  gt_queue_delete(cuss->outqueue);
  for (i = 0; i < gt_array_size(cuss->cur_cds_set); i++) {
    gt_genome_node_delete(*(GtGenomeNode**) gt_array_get(cuss->cur_cds_set, i));
  }
  gt_node_stream_delete(cuss->in_stream);
  gt_array_delete(cuss->cur_cds_set);
}

const GtNodeStreamClass* gt_codon_usage_scan_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtCodonUsageScanStream),
                                   codon_usage_scan_stream_free,
                                   codon_usage_scan_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

GtNodeStream* gt_codon_usage_scan_stream_new(GtNodeStream *in_stream)
{
  GtCodonUsageScanStream *cuss;
  GtNodeStream *ns;
  gt_assert(in_stream);
  ns = gt_node_stream_create(gt_codon_usage_scan_stream_class(), true);
  cuss = gt_codon_usage_scan_stream_cast(ns);
  cuss->cur_cds_set = gt_array_new(sizeof (GtFeatureNode*));
  cuss->in_stream = gt_node_stream_ref(in_stream);
  cuss->cur_cds_range.start =
    cuss->cur_cds_range.end = GT_UNDEF_UWORD;
  cuss->outqueue = gt_queue_new();
  return ns;
}
