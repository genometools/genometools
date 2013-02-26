/*
  Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "core/array.h"
#include "core/class_alloc_lock.h"
#include "core/cstr_api.h"
#include "core/encseq_api.h"
#include "core/hashmap.h"
#include "core/log.h"
#include "core/str_api.h"
#include "core/str_array.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/feature_type.h"
#include "extended/node_visitor_api.h"
#include "ltr/ltr_cluster_prepare_seq_visitor.h"

struct GtLTRClusterPrepareSeqVisitor {
  const GtNodeVisitor parent_instance;
  GtEncseq *src_encseq;
  GtHashmap *feat_to_encseq,
            *encseq_builders;
  GtStrArray *feat_to_encseq_keys;
};

static int extract_feature_seq(GtEncseqBuilder *b, const char *header,
                               GtStr *seqid, GtEncseq *encseq, GtRange range,
                               GT_UNUSED const char *fnt, GtError *err)
{
  char *buffer;
  int had_err = 0;
  unsigned long seqnum,
                startpos;

  (void) sscanf(gt_str_get(seqid), "seq%lu", &seqnum);
  if (seqnum >= gt_encseq_num_of_sequences(encseq)) {
    gt_error_set(err, "annotation encountered for sequence %lu, but the "
                      "supplied encoded sequence only contains sequences "
                      "0-%lu", seqnum, gt_encseq_num_of_sequences(encseq)-1);
    had_err = -1;
  }

  if (!had_err) {
    buffer = gt_calloc((size_t) gt_range_length(&range) + 1, sizeof (char));
    startpos = gt_encseq_seqstartpos(encseq, seqnum);
    gt_encseq_extract_decoded(encseq, buffer, startpos + range.start,
                              startpos + range.end);
    gt_encseq_builder_add_cstr(b, buffer, gt_range_length(&range), header);
    gt_free(buffer);
  }

  return had_err;
}

static int gt_ltr_cluster_prepare_seq_visitor_feature_node(GtNodeVisitor *nv,
                                               GtFeatureNode *fn,
                                               GtError *err)
{
  GtLTRClusterPrepareSeqVisitor *lcv;
  GtFeatureNode *curnode;
  GtFeatureNodeIterator *fni;
  GtStr *seqid = NULL;
  const char *fnt;
  char buffer[BUFSIZ];
  int had_err = 0;
  bool first_ltr = true;
  gt_error_check(err);

  lcv = gt_ltr_cluster_prepare_seq_visitor_cast(nv);
  gt_assert(lcv);
  gt_error_check(err);
  fni = gt_feature_node_iterator_new(fn);

  while (!had_err && (curnode = gt_feature_node_iterator_next(fni))) {
    GtFile *file = NULL;
    fnt = gt_feature_node_get_type(curnode);
    if (strcmp(fnt, gt_ft_repeat_region) == 0) {
      const char *rid;
      unsigned long id;
      seqid = gt_genome_node_get_seqid((GtGenomeNode*) curnode);
      rid = gt_feature_node_get_attribute(curnode, "ID");
      (void) sscanf(rid, "repeat_region%lu", &id);
      (void) snprintf(buffer, BUFSIZ, "%s_%lu", gt_str_get(seqid), id);
    } else if (strcmp(fnt, gt_ft_protein_match) == 0) {
      GtRange range;
      const char *attr;
      GtEncseqBuilder *eb;
      char header[BUFSIZ];
      attr = gt_feature_node_get_attribute(curnode, "name");
      if (!attr)
        continue;
      range = gt_genome_node_get_range((GtGenomeNode*) curnode);
      (void) snprintf(header, BUFSIZ, "%s_%lu_%lu", buffer, range.start,
                      range.end);
      if (!gt_hashmap_get(lcv->encseq_builders, attr)) {
        eb = gt_encseq_builder_new(gt_encseq_alphabet(lcv->src_encseq));
        gt_encseq_builder_create_ssp_tab(eb);
        gt_encseq_builder_create_sds_tab(eb);
        gt_encseq_builder_create_des_tab(eb);
        gt_hashmap_add(lcv->encseq_builders, gt_cstr_dup(attr),
                       eb);
        gt_log_log("builder %p added for feature %s", eb, attr);
        gt_str_array_add_cstr(lcv->feat_to_encseq_keys, attr);
      } else {
        eb = (GtEncseqBuilder*) gt_hashmap_get(lcv->encseq_builders, attr);
      }
      had_err = extract_feature_seq(eb, header, seqid, lcv->src_encseq, range,
                                    fnt, err);
      gt_file_delete(file);
    } else if (strcmp(fnt, gt_ft_LTR_retrotransposon) == 0)
      continue;
    else {
      char *tmp;
      GtRange range;
      GtEncseqBuilder *eb;
      char header[BUFSIZ];  /* XXX: use GtStr for safety */
      if (strcmp(fnt, gt_ft_long_terminal_repeat) == 0) {
        if (first_ltr) {
          tmp = gt_cstr_dup("lLTR");
          first_ltr = false;
        } else
          tmp = gt_cstr_dup("rLTR");
      } else
        tmp = gt_cstr_dup(fnt);
      range = gt_genome_node_get_range((GtGenomeNode*) curnode);
      if ((range.end - range.start + 1) < 10UL) {
        gt_free(tmp);
        continue;
      }
      (void) snprintf(header, BUFSIZ, "%s_%lu_%lu", buffer, range.start,
                      range.end);
      if (!gt_hashmap_get(lcv->encseq_builders, tmp)) {
        eb = gt_encseq_builder_new(gt_encseq_alphabet(lcv->src_encseq));
        gt_encseq_builder_create_ssp_tab(eb);
        gt_encseq_builder_create_sds_tab(eb);
        gt_encseq_builder_create_des_tab(eb);
        gt_hashmap_add(lcv->encseq_builders, gt_cstr_dup(tmp),
                       eb);
        gt_str_array_add_cstr(lcv->feat_to_encseq_keys, tmp);
      } else {
        eb = (GtEncseqBuilder*) gt_hashmap_get(lcv->encseq_builders, tmp);
      }
      had_err = extract_feature_seq(eb, header, seqid, lcv->src_encseq, range,
                                    fnt, err);
      gt_free(tmp);
      gt_file_delete(file);
    }
  }
  gt_feature_node_iterator_delete(fni);

  return had_err;
}

void gt_ltr_cluster_prepare_seq_visitor_free(GtNodeVisitor *v)
{
  GtLTRClusterPrepareSeqVisitor *lcv;
  lcv = gt_ltr_cluster_prepare_seq_visitor_cast(v);
  gt_encseq_delete(lcv->src_encseq);
  gt_str_array_delete(lcv->feat_to_encseq_keys);
  gt_hashmap_delete(lcv->feat_to_encseq);
  gt_hashmap_delete(lcv->encseq_builders);
}

const GtNodeVisitorClass* gt_ltr_cluster_prepare_seq_visitor_class(void)
{
  static const GtNodeVisitorClass *nvc;
  gt_class_alloc_lock_enter();
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof(GtLTRClusterPrepareSeqVisitor),
                                gt_ltr_cluster_prepare_seq_visitor_free,
                                NULL,
                                gt_ltr_cluster_prepare_seq_visitor_feature_node,
                                NULL,
                                NULL,
                                NULL);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}

GtNodeVisitor* gt_ltr_cluster_prepare_seq_visitor_new(GtEncseq *encseq,
                                                      GT_UNUSED GtError *err)
{
  GtNodeVisitor *nv;
  GtLTRClusterPrepareSeqVisitor *lcv;
  nv = gt_node_visitor_create(gt_ltr_cluster_prepare_seq_visitor_class());
  lcv = gt_ltr_cluster_prepare_seq_visitor_cast(nv);
  gt_assert(lcv);
  lcv->src_encseq = gt_encseq_ref(encseq);
  lcv->feat_to_encseq = NULL;
  lcv->encseq_builders = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                        (GtFree) gt_encseq_builder_delete);
  lcv->feat_to_encseq_keys = gt_str_array_new();
  return nv;
}

int gt_ltr_cluster_prepare_seq_finish_encseqs(void *key, void *value,
                                              void *data, GtError *err)
{
  GtLTRClusterPrepareSeqVisitor *v = (GtLTRClusterPrepareSeqVisitor*) data;
  GtEncseqBuilder *eb = (GtEncseqBuilder*) value;
  GtEncseq *encseq;
  int had_err = 0;
  const char *feature = (const char*) key;
  gt_assert(key && value);
  gt_error_check(err);

  encseq = gt_encseq_builder_build(eb, err);
  if (!encseq)
    had_err = -1;
  if (!had_err) {
    if (!v->feat_to_encseq) {
      v->feat_to_encseq = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                         (GtFree) gt_encseq_delete);
    }
    gt_assert(!gt_hashmap_get(v->feat_to_encseq, feature));
    gt_hashmap_add(v->feat_to_encseq, gt_cstr_dup(feature), encseq);
    gt_log_log("added encseq %p to hash for feature %s", encseq, feature);
  }
  return had_err;
}

GtHashmap* gt_ltr_cluster_prepare_seq_visitor_get_encseqs(
                                              GtLTRClusterPrepareSeqVisitor *v)
{
  gt_assert(v && v->encseq_builders);
  gt_log_log("finishing encseqs");
  if (!v->feat_to_encseq) {
    gt_log_log("starting...");
    (void) gt_hashmap_foreach(v->encseq_builders,
                              gt_ltr_cluster_prepare_seq_finish_encseqs,
                              v,
                              NULL);
  }
  return gt_hashmap_ref(v->feat_to_encseq);
}

GtStrArray* gt_ltr_cluster_prepare_seq_visitor_get_features(
                                              GtLTRClusterPrepareSeqVisitor *v)
{
  gt_assert(v && v->feat_to_encseq_keys);
  return gt_str_array_ref(v->feat_to_encseq_keys);
}
