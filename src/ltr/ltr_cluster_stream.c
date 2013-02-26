/*
  Copyright (c) 2011      Sascha Kastens <sascha.kastens@studium.uni-hamburg.de>
  Copyright (c)      2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2011-2012 Center for Bioinformatics, University of Hamburg

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

#include "core/array.h"
#include "core/class_alloc_lock.h"
#include "core/cstr_api.h"
#include "core/encseq.h"
#include "core/hashmap.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/str_array.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/clustered_set.h"
#include "extended/clustered_set_uf.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/feature_type.h"
#include "extended/node_stream_api.h"
#include "extended/match.h"
#include "extended/match_iterator_api.h"
#include "extended/match_iterator_last.h"
#include "extended/match_iterator_open.h"
#include "ltr/ltr_cluster_stream.h"
#include "ltr/ltr_cluster_prepare_seq_visitor.h"
#include "match/sfx-run.h"

struct GtLTRClusterStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtLTRClusterPrepareSeqVisitor *lcv;
  GtArray *nodes;
  GtHashmap *feat_to_encseq;
  GtStrArray *feat_to_encseq_keys;
  bool first_next;
  unsigned long psmall,
                plarge,
                next_index;
  int match_score, mismatch_cost, gap_open_cost,
      gap_ext_cost, xdrop, ydrop, zdrop, mscoregapped,
      mscoregapless, k;
  char **current_state;
};

#define gt_ltr_cluster_stream_cast(CS)\
        gt_node_stream_cast(gt_ltr_cluster_stream_class(), CS);

typedef struct GtMatchReference {
  unsigned long startpos, matchnum;
} GtMatchReference;

typedef struct GtMatchEdge {
unsigned long matchnum0, matchnum1;
unsigned long gap_size;
unsigned long edist;
unsigned long minlength;
} GtMatchEdge;

typedef struct GtMatchEdgeTable {
  GtArray *edges;
  unsigned long num_of_edges;
} GtMatchEdgeTable;

#define STORECLUSTEREDGEG(I, J, VALUE)\
matchedge.matchnum0 = I;\
matchedge.matchnum1 = J;\
matchedge.gap_size = VALUE;\
gt_array_add(matchedgetab.edges, matchedge);\
matchedgetab.num_of_edges++

#define STORECLUSTEREDGEED(I, J, ML, ED)\
matchedge.matchnum0 = I;\
matchedge.matchnum1 = J;\
matchedge.minlength = ML;\
matchedge.edist = ED;\
gt_array_add(matchedgetab.edges, matchedge);\
matchedgetab.num_of_edges++

static int cmpmatchreferences(const void *dataA, const void *dataB)
{
  GtMatchReference *mrefA = (GtMatchReference*) dataA;
  GtMatchReference *mrefB = (GtMatchReference*) dataB;
  if (mrefA->startpos < mrefB->startpos)
    return -1;
  if (mrefA->startpos > mrefB->startpos)
    return 1;
  return 0;
}

static GtMatchReference* gt_mirror_and_sort_matches(GtArray *matches)
{
  GtMatchReference *mref;
  GtMatch *match;
  GtRange rng_seq1, rng_seq2;
  unsigned long i, j;
  mref = gt_calloc((size_t) (gt_array_size(matches) * 2),
                   sizeof (GtMatchReference));
  for (i = 0, j = 0; i < gt_array_size(matches); i++, j+=2) {
    match = *(GtMatch**) gt_array_get(matches, i);
    gt_match_get_range_seq1(match, &rng_seq1);
    gt_match_get_range_seq2(match, &rng_seq2);
    mref[j].startpos = rng_seq1.start;
    mref[j].matchnum = i;
    mref[j+1].startpos = rng_seq2.start;
    mref[j+1].matchnum = i;
  }
  qsort (mref,
         (size_t) (2 * gt_array_size(matches)),
         sizeof (GtMatchReference),
         cmpmatchreferences);
  return mref;
}

static int gt_cluster_matches(GtClusteredSet *cs,
                              GtMatchEdgeTable *matchedgetab,
                              GtError *err)
{
  unsigned long i;
  for (i = 0; i < matchedgetab->num_of_edges; i++) {
     if (gt_clustered_set_merge_clusters(cs,
              ((GtMatchEdge*)(gt_array_get(matchedgetab->edges, i)))->matchnum0,
              ((GtMatchEdge*)(gt_array_get(matchedgetab->edges, i)))->matchnum1,
              err) != 0) {
       return -1;
    }
  }
  return 0;
}

static int cluster_sequences(GtArray *matches,
                             GtClusteredSet *cs,
                             GtHashmap *seqdesc2seqnum,
                             unsigned int psmall,
                             unsigned int plarge,
                             GtEncseq *encseq,
                             GtError *err)
{
  GtMatch *match;
  GtMatchEdgeTable matchedgetab;
  GtMatchEdge matchedge;
  GtRange rng_seq1,
          rng_seq2;
  int had_err = 0;
  unsigned long i,
                lsmall,
                llarge,
                matchlen1,
                matchlen2,
                num_of_seq,
                seqnum1 = 0,
                seqnum2 = 0;
  const char *seqid;

  num_of_seq = gt_encseq_num_of_sequences(encseq);
  gt_assert(matches && cs && seqdesc2seqnum && encseq);

  if (gt_clustered_set_num_of_elements(cs, err) != num_of_seq) {
    had_err = -1;
    gt_error_set(err,
                 "number of sequences (%lu) unequals number of elements in"
                 " clustered set (%lu)",
                 num_of_seq, gt_clustered_set_num_of_elements(cs, err));
  }
  if (!had_err) {
    matchedgetab.edges = gt_array_new(sizeof (GtMatchEdge));
    matchedgetab.num_of_edges = 0;

    for (i = 0; i < gt_array_size(matches); i++) {
      match = *(GtMatch**) gt_array_get(matches, i);
      gt_match_get_range_seq1(match, &rng_seq1);
      gt_match_get_range_seq2(match, &rng_seq2);

      matchlen1 =  gt_range_length(&rng_seq1);
      matchlen2 =  gt_range_length(&rng_seq2);

      seqid = gt_match_get_seqid1(match);
      if (gt_hashmap_get(seqdesc2seqnum, (void*) seqid) != NULL)
        seqnum1 = ((unsigned long) gt_hashmap_get(seqdesc2seqnum, seqid)) - 1;
      else {
        had_err = -1;
        gt_error_set(err, "key %s not found", seqid);
      }

      seqid = gt_match_get_seqid2(match);
      if (!had_err && gt_hashmap_get(seqdesc2seqnum, (void*) seqid) != NULL)
        seqnum2 = ((unsigned long) gt_hashmap_get(seqdesc2seqnum, seqid)) - 1;
      else {
        had_err = -1;
        gt_error_set(err, "key %s not found", seqid);
      }

      if (!had_err) {
        if (gt_encseq_seqlength(encseq, seqnum1) >
            gt_encseq_seqlength(encseq, seqnum2)) {
          llarge = gt_encseq_seqlength(encseq, seqnum1);
          lsmall = gt_encseq_seqlength(encseq, seqnum2);
        } else {
          lsmall = gt_encseq_seqlength(encseq, seqnum1);
          llarge = gt_encseq_seqlength(encseq, seqnum2);
        }
        if (((llarge * plarge)/100 <= matchlen1) &&
            ((lsmall * psmall)/100 <= matchlen1) &&
            ((llarge * plarge)/100 <= matchlen2) &&
            ((lsmall * psmall)/100 <= matchlen2)) {
          if (seqnum1 != seqnum2) {
            matchedge.matchnum0 = seqnum1;
            matchedge.matchnum1 = seqnum2;
            gt_array_add(matchedgetab.edges, matchedge);
            matchedgetab.num_of_edges++;
          }
        }
      }
    }
  }
  if (!had_err)
    if (gt_cluster_matches(cs, &matchedgetab, err) != 0)
      had_err = -1;
  if (!had_err)
    gt_array_delete(matchedgetab.edges);
  return had_err;
}

GT_UNUSED static int gt_cluster_matches_gap(GtArray *matches,
                                            GtClusteredSet *cs,
                                            unsigned long max_gap_size,
                                            GtError *err)
{
  GtMatchReference *mref;
  GtMatchEdgeTable matchedgetab;
  GtMatchEdge matchedge;
  GtMatch *match;
  GtRange range;
  unsigned long i,
                j,
                gap_size = 0,
                length = 0;
  unsigned long num_of_matches;
  int had_err;

  num_of_matches = gt_array_size(matches);
  if (gt_clustered_set_num_of_elements(cs, err) != num_of_matches) {
    had_err = -1;
    gt_error_set(err,
                 "number of matches (%lu) unequals number of elements in"
                 "clustered set (%lu)",
                 num_of_matches, gt_clustered_set_num_of_elements(cs, err));
  }
  if (!had_err) {
    matchedgetab.edges = gt_array_new(sizeof (GtMatchEdge));
    matchedgetab.num_of_edges = 0;

    mref = gt_mirror_and_sort_matches(matches);

    for (i = 0; i < (2 * (num_of_matches) - 1); i++) {
      for (j = i + 1; j < (2 * num_of_matches); j++) {
        match = *(GtMatch**) gt_array_get(matches, mref[i].matchnum);
        gt_match_get_range_seq1(match, &range);
        length = gt_range_length(&range);
        gap_size = mref[j].startpos - mref[i].startpos + length;

        if (gap_size > max_gap_size)
          break;
        if (mref[i].matchnum != mref[j].matchnum) {
          STORECLUSTEREDGEG(mref[i].matchnum, mref[j].matchnum, gap_size);
        }
      }
    }
    if (gt_cluster_matches(cs, &matchedgetab, err) != 0)
      had_err = -1;
    gt_array_delete(matchedgetab.edges);
    gt_free(mref);
  }
  return had_err;
}

static void free_hash(void *elem)
{
  gt_free(elem);
}

static int cluster_annotate_nodes(GtClusteredSet *cs, GtEncseq *encseq,
                                  const char *feature, GtArray *nodes,
                                  GtError *err)
{
  GtFeatureNodeIterator *fni;
  GtFeatureNode *curnode = NULL, *tmp;
  GtClusteredSetIterator *csi = NULL;
  GtGenomeNode *gn;
  GtHashmap *desc2node;
  GtStr *seqid = NULL;
  int had_err = 0;
  unsigned long num_of_clusters, i, elm;
  const char *fnt = NULL;
  char buffer[BUFSIZ], *real_feature;
  gt_error_check(err);

  if ((strcmp(feature, "lLTR") == 0) || (strcmp(feature, "rLTR") == 0))
    real_feature = gt_cstr_dup(gt_ft_long_terminal_repeat);
  else
    real_feature = gt_cstr_dup(feature);

  desc2node = gt_hashmap_new(GT_HASH_STRING, free_hash, NULL);
  for (i = 0; i < gt_array_size(nodes); i++) {
    gn = *(GtGenomeNode**) gt_array_get(nodes, i);
    if (gt_feature_node_try_cast(gn) == NULL)
      continue;
    fni = gt_feature_node_iterator_new((GtFeatureNode*) gn);
    while ((curnode = gt_feature_node_iterator_next(fni)) != NULL) {
      char header[BUFSIZ];
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
        attr = gt_feature_node_get_attribute(curnode, "name");
        if (!attr)
          continue;
        if (strcmp(feature, attr) != 0)
          continue;
        range = gt_genome_node_get_range((GtGenomeNode*) curnode);
        if ((range.end - range.start + 1) < 10UL)
          continue;
        (void) snprintf(header, BUFSIZ, "%s_%lu_%lu", buffer, range.start,
                        range.end);
        gt_hashmap_add(desc2node, (void*) gt_cstr_dup(header), (void*) curnode);
      } else if (strcmp(fnt, real_feature) == 0) {
        GtRange range;
        range = gt_genome_node_get_range((GtGenomeNode*) curnode);
        if ((range.end - range.start + 1) < 10UL)
          continue;
        (void) snprintf(header, BUFSIZ, "%s_%lu_%lu", buffer, range.start,
                        range.end);
        gt_hashmap_add(desc2node, (void*) gt_cstr_dup(header), (void*) curnode);
      }
    }
    gt_feature_node_iterator_delete(fni);
  }
  gt_free(real_feature);

  num_of_clusters = gt_clustered_set_num_of_clusters(cs, err);
  for (i = 0; i < num_of_clusters; i++) {
    csi = gt_clustered_set_get_iterator(cs, i ,err);
    if (csi != NULL) {
      while (!had_err && (gt_clustered_set_iterator_next(csi, &elm, err)
             != GT_CLUSTERED_SET_ITERATOR_STATUS_END)) {
        char clid[BUFSIZ];
        const char *encseqdesc;
        char *encseqid;
        unsigned long desclen;
        encseqdesc = gt_encseq_description(encseq, &desclen, elm);
        encseqid = gt_calloc((size_t) (desclen + 1), sizeof (char));
        (void) strncpy(encseqid, encseqdesc, (size_t) desclen);
        encseqid[desclen] = '\0';
        tmp = (GtFeatureNode*) gt_hashmap_get(desc2node, (void*) encseqid);
        (void) snprintf(clid, BUFSIZ, "%lu", i);
        gt_feature_node_set_attribute(tmp, "clid", clid);
        gt_free(encseqid);
      }
    }
    gt_clustered_set_iterator_delete(csi, err);
    csi = NULL;
  }
  gt_hashmap_delete(desc2node);
  return had_err;
}

static int process_feature(GtLTRClusterStream *lcs,
                           const char *feature,
                           GtError *err)
{
  GtArray *matches;
  GtMatchIterator *mi = NULL;
  GtMatch *match = NULL;
  GtMatchIteratorStatus status;
  GtEncseq *encseq;
  unsigned long i;
  int had_err = 0;

  if (lcs->current_state != NULL) {
    char tmp[BUFSIZ];
    gt_free(*lcs->current_state);
    (void) snprintf(tmp, BUFSIZ, "Clustering feature: %s", feature);
    *lcs->current_state = gt_cstr_dup(tmp);
  }
  matches = gt_array_new(sizeof(GtMatch*));
  encseq = (GtEncseq*) gt_hashmap_get(lcs->feat_to_encseq, feature);
  gt_log_log("found encseq %p for feature %s", encseq, feature);
  if (!had_err) {
    mi = gt_match_iterator_last_new(encseq, encseq, lcs->match_score,
                                        lcs->mismatch_cost,
                                        lcs->gap_open_cost,
                                        lcs->gap_ext_cost,
                                        lcs->xdrop,
                                        lcs->ydrop,
                                        lcs->zdrop,
                                        lcs->k,
                                        lcs->mscoregapped,
                                        lcs->mscoregapless, err);
    if (mi != NULL) {
      while ((status = gt_match_iterator_next(mi, &match, err))
             != GT_MATCHER_STATUS_END) {
        if (status == GT_MATCHER_STATUS_OK) {
          gt_array_add(matches, match);
        } else {
          gt_assert(status == GT_MATCHER_STATUS_ERROR);
          had_err = -1;
          break;
        }
      }
    } else
      had_err = -1;
  }
  if (!had_err) {
    GtClusteredSet *cs;
    GtHashmap *seqdesc2seqnum;
    GtMatch *tmp_match;
    const char *description;
    char *output;
    unsigned long desclen,
                  num_of_seq;

    seqdesc2seqnum = gt_hashmap_new(GT_HASH_STRING, free_hash, NULL);
    num_of_seq = gt_encseq_num_of_sequences(encseq);
    for (i = 0; i < num_of_seq; i++) {
      description = gt_encseq_description(encseq, &desclen, i);
      output = gt_calloc((size_t) (desclen + 1), sizeof (char));
      strncpy(output, description, (size_t) desclen);
      output[desclen] = '\0';
      gt_hashmap_add(seqdesc2seqnum, (void*) gt_cstr_dup(output),
                     (void*) (i + 1));
      gt_free(output);
    }
    cs = gt_clustered_set_union_find_new(num_of_seq, err);
    if (cs != NULL) {
      if (cluster_sequences(matches, cs, seqdesc2seqnum, (unsigned) lcs->psmall,
                            (unsigned) lcs->plarge, encseq, err) != 0) {
        had_err = -1;
      }
      if (!had_err) {
        (void) cluster_annotate_nodes(cs, encseq, feature, lcs->nodes, err);
      }
    } else
      had_err = -1;

    for (i = 0; i < gt_array_size(matches); i++) {
      tmp_match = *(GtMatch**) gt_array_get(matches, i);
      gt_match_delete(tmp_match);
    }
    gt_array_delete(matches);
    matches = NULL;
    gt_hashmap_delete(seqdesc2seqnum);
    gt_clustered_set_delete(cs, err);
  }
  gt_match_iterator_delete(mi);

  return had_err;
}

static int gt_ltr_cluster_stream_next(GtNodeStream *ns,
                                      GtGenomeNode **gn,
                                      GtError *err)
{
  GtLTRClusterStream *lcs;
  GtGenomeNode *ref_gn;
  int had_err = 0;
  unsigned long i = 0;

  gt_error_check(err);
  lcs = gt_ltr_cluster_stream_cast(ns);
  if (lcs->first_next) {
    while (!(had_err = gt_node_stream_next(lcs->in_stream, gn, err)) && *gn) {
      gt_assert(*gn && !had_err);
      ref_gn = gt_genome_node_ref(*gn);
      gt_array_add(lcs->nodes, ref_gn);
      had_err = gt_genome_node_accept(*gn, (GtNodeVisitor*) lcs->lcv, err);
      if (had_err) {
        gt_genome_node_delete(*gn);
        *gn = NULL;
        break;
      }
    }
    lcs->feat_to_encseq =
                       gt_ltr_cluster_prepare_seq_visitor_get_encseqs(lcs->lcv);
    lcs->feat_to_encseq_keys =
                      gt_ltr_cluster_prepare_seq_visitor_get_features(lcs->lcv);
    if (!had_err) {
      for (i = 0; i < gt_str_array_size(lcs->feat_to_encseq_keys); i++) {
        had_err = process_feature(lcs,
                                  gt_str_array_get(lcs->feat_to_encseq_keys, i),
                                  err);
        if (had_err)
          break;
      }
    }
    if (!had_err) {
      *gn = *(GtGenomeNode**) gt_array_get(lcs->nodes, lcs->next_index);
      lcs->next_index++;
      lcs->first_next = false;
      return 0;
    }
  } else {
    if (lcs->next_index >= gt_array_size(lcs->nodes))
      *gn = NULL;
    else {
      *gn = *(GtGenomeNode**) gt_array_get(lcs->nodes, lcs->next_index);
      lcs->next_index++;
    }
    return 0;
  }

  return had_err;
}

static void gt_ltr_cluster_stream_free(GtNodeStream *ns)
{
  unsigned long i;
  GtLTRClusterStream *lcs = gt_ltr_cluster_stream_cast(ns);
  gt_node_visitor_delete((GtNodeVisitor*) lcs->lcv);
  gt_hashmap_delete(lcs->feat_to_encseq);
  gt_str_array_delete(lcs->feat_to_encseq_keys);
  for (i = 0; i < gt_array_size(lcs->nodes); i++)
    gt_genome_node_delete(*(GtGenomeNode**) gt_array_get(lcs->nodes, i));
  gt_array_delete(lcs->nodes);
  gt_node_stream_delete(lcs->in_stream);
}

const GtNodeStreamClass* gt_ltr_cluster_stream_class(void)
{
  static const GtNodeStreamClass *nsc;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof(GtLTRClusterStream),
                                   gt_ltr_cluster_stream_free,
                                   gt_ltr_cluster_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

GtNodeStream* gt_ltr_cluster_stream_new(GtNodeStream *in_stream,
                                        GtEncseq *encseq,
                                        int match_score,
                                        int mismatch_cost,
                                        int gap_open_cost,
                                        int gap_ext_cost,
                                        int xdrop,
                                        int ydrop,
                                        int zdrop,
                                        int k,
                                        int mscoregapped,
                                        int mscoregapless,
                                        unsigned long plarge,
                                        unsigned long psmall,
                                        char **current_state,
                                        GtError *err)
{
  GtNodeStream *ns;
  GtLTRClusterStream *lcs;
  ns = gt_node_stream_create(gt_ltr_cluster_stream_class(), false);
  lcs = gt_ltr_cluster_stream_cast(ns);
  lcs->in_stream = gt_node_stream_ref(in_stream);
  lcs->feat_to_encseq = NULL;
  lcs->feat_to_encseq_keys = NULL;
  lcs->nodes = gt_array_new(sizeof(GtGenomeNode*));
  lcs->lcv = gt_ltr_cluster_prepare_seq_visitor_cast(
                           gt_ltr_cluster_prepare_seq_visitor_new(encseq, err));
  lcs->first_next = true;
  lcs->next_index = 0;
  lcs->match_score = match_score;
  lcs->mismatch_cost = mismatch_cost;
  lcs->gap_open_cost = gap_open_cost;
  lcs->gap_ext_cost = gap_ext_cost;
  lcs->xdrop = xdrop;
  lcs->ydrop = ydrop;
  lcs->zdrop = zdrop;
  lcs->mscoregapped = mscoregapped;
  lcs->mscoregapless = mscoregapless;
  lcs->k = k;
  lcs->plarge = plarge;
  lcs->psmall = psmall;
  lcs->current_state = current_state;
  return ns;
}
