/*
  Copyright (c) 2011 Sascha Kastens <sascha.kastens@studium.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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
#include "core/cstr_api.h"
#include "core/encseq_api.h"
#include "core/fasta.h"
#include "core/hashmap.h"
#include "core/str_api.h"
#include "core/str_array.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/node_visitor_api.h"
#include "ltr/ltr_cluster_visitor.h"

struct GtLTRClusterVisitor {
  const GtNodeVisitor parent_instance;
  GtEncseq *encseq;
  GtHashmap *feat_to_file;
  GtStrArray *feat_to_file_keys;
  GtStr *file_prefix;
};

static int gt_ltr_cluster_visitor_extract_feature_seq(GtFile *file,
                               const char *header, GtStr *seqid,
                               GtEncseq *encseq, GtRange range)
{
  char *buffer;
  int had_err = 0;
  unsigned long seqnum,
                startpos;

  (void) sscanf(gt_str_get(seqid), "seq%lu", &seqnum);
  buffer = gt_calloc((size_t) gt_range_length(&range) + 1, sizeof (char));
  startpos = gt_encseq_seqstartpos(encseq, seqnum);
  gt_encseq_extract_decoded(encseq, buffer, startpos + range.start - 1,
                            startpos + range.end - 1);
  gt_fasta_show_entry(header, buffer, gt_range_length(&range), 50UL, file);
  gt_free(buffer);

  return had_err;

}

static int gt_ltr_cluster_visitor_feature_node(GtNodeVisitor *nv,
                                               GtFeatureNode *fn,
                                               GtError *err)
{
  GtLTRClusterVisitor *lcv;
  GtFeatureNode *curnode;
  GtFeatureNodeIterator *fni;
  GtStr *seqid = NULL;
  const char *fnt;
  char buffer[BUFSIZ];
  int had_err = 0;
  bool first_ltr = true;

  lcv = gt_ltr_cluster_visitor_cast(nv);
  gt_assert(lcv);
  gt_error_check(err);
  fni = gt_feature_node_iterator_new(fn);

  while (!had_err && (curnode = gt_feature_node_iterator_next(fni))) {
    GtFile *file = NULL;
    fnt = gt_feature_node_get_type(curnode);
    if (strcmp(fnt, "repeat_region") == 0) {
      const char *rid;
      unsigned long id;
      seqid = gt_genome_node_get_seqid((GtGenomeNode*) curnode);
      rid = gt_feature_node_get_attribute(curnode, "ID");
      (void) sscanf(rid, "repeat_region%lu", &id);
      (void) snprintf(buffer, BUFSIZ, "%s_%lu", gt_str_get(seqid), id);
    } else if (strcmp(fnt, "protein_match") == 0) {
      GtRange range;
      const char *attr;
      char header[BUFSIZ];
      attr = gt_feature_node_get_attribute(curnode, "name");
      if (!attr)
        continue;
      range = gt_genome_node_get_range((GtGenomeNode*) curnode);
      (void) snprintf(header, BUFSIZ, "%s_%lu_%lu", buffer, range.start,
                      range.end);
      if (!gt_hashmap_get(lcv->feat_to_file, (void*) attr)) {
        char filename[BUFSIZ];
        (void) snprintf(filename, BUFSIZ, "%s_%s.fas",
                        gt_str_get(lcv->file_prefix), attr);
        gt_hashmap_add(lcv->feat_to_file, (void*) gt_cstr_dup(attr),
                       (void*) gt_cstr_dup(filename));
        gt_str_array_add_cstr(lcv->feat_to_file_keys, attr);
        file = gt_file_new(filename, "w", err);
      } else {
        char *fname;
        fname = (char*) gt_hashmap_get(lcv->feat_to_file, (void*) attr);
        file = gt_file_new(fname, "a", err);
      }
      had_err = gt_ltr_cluster_visitor_extract_feature_seq(file, header,
                                                 seqid, lcv->encseq, range);
      gt_file_delete(file);
    } else if (strcmp(fnt, "LTR_retrotransposon") == 0)
      continue;
    else {
      char *tmp;
      GtRange range;
      char header[BUFSIZ];
      if (strcmp(fnt, "long_terminal_repeat") == 0) {
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
      if (!gt_hashmap_get(lcv->feat_to_file, (void*) tmp)) {
        char filename[BUFSIZ];
        (void) snprintf(filename, BUFSIZ, "%s_%s.fas",
                        gt_str_get(lcv->file_prefix), tmp);
        gt_hashmap_add(lcv->feat_to_file, (void*) gt_cstr_dup(tmp),
                       (void*) gt_cstr_dup(filename));
        gt_str_array_add_cstr(lcv->feat_to_file_keys, tmp);
        file = gt_file_new(filename, "w", err);
      } else {
        char *fname;
        fname = (char*) gt_hashmap_get(lcv->feat_to_file, (void*) tmp);
        file = gt_file_new(fname, "a", err);
      }
      had_err = gt_ltr_cluster_visitor_extract_feature_seq(file, header,
                                                    seqid, lcv->encseq, range);
      gt_free(tmp);
      gt_file_delete(file);
    }
  }
  gt_feature_node_iterator_delete(fni);

  return had_err;
}

const GtNodeVisitorClass* gt_ltr_cluster_visitor_class(void)
{
  static const GtNodeVisitorClass *nvc;
  if (!nvc)
    nvc = gt_node_visitor_class_new(sizeof(GtLTRClusterVisitor),
                                    NULL,
                                    NULL,
                                    gt_ltr_cluster_visitor_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  return nvc;
}

GtNodeVisitor* gt_ltr_cluster_visitor_new(GtEncseq *encseq,
                                          GtStr *file_prefix,
                                          GtHashmap *feat_to_file,
                                          GtStrArray *feat_to_file_keys,
                                          GT_UNUSED GtError *err)
{
  GtNodeVisitor *nv;
  GtLTRClusterVisitor *lcv;
  nv = gt_node_visitor_create(gt_ltr_cluster_visitor_class());
  lcv = gt_ltr_cluster_visitor_cast(nv);
  gt_assert(lcv);
  lcv->encseq = encseq;
  lcv->feat_to_file = feat_to_file;
  lcv->feat_to_file_keys = feat_to_file_keys;
  lcv->file_prefix = file_prefix;
  return nv;
}
