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

#include "core/array.h"
#include "core/cstr_api.h"
#include "core/encseq.h"
#include "core/hashmap.h"
#include "core/ma.h"
#include "core/str_array.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/clustered_set.h"
#include "extended/clustered_set_uf.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/node_stream_api.h"
#include "extended/match.h"
#include "extended/match_iterator_api.h"
#include "extended/match_iterator_blast.h"
#include "extended/match_iterator_open.h"
#include "ltr/ltr_cluster_stream.h"
#include "ltr/ltr_cluster_visitor.h"
#include "match/sfx-run.h"

struct GtLTRClusterStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtLTRClusterVisitor *lcv;
  GtArray *nodes;
  GtHashmap *feat_to_file;
  GtStrArray *feat_to_file_keys;
  bool dust;
  int word_size,
      gapopen,
      gapextend,
      penalty,
      reward,
      num_threads;
  double evalue,
         identity,
         xdrop;
  bool first_next,
       from_file;
  unsigned long psmall,
                plarge,
                next_index;
  const char *moreblast;
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
          STORECLUSTEREDGEG (mref[i].matchnum, mref[j].matchnum, gap_size);
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
    real_feature = gt_cstr_dup("long_terminal_repeat");
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

static void write_match_to_file(GtMatch *match, GtFile *matchf)
{
  char buffer[BUFSIZ];
  GtRange rng1,
          rng2;

  gt_match_get_range_seq1(match, &rng1);
  gt_match_get_range_seq2(match, &rng2);
  (void) snprintf(buffer, BUFSIZ, " %lu %s %lu c %lu %s %lu 0 0.0e0 0 0.0\n",
                  gt_range_length(&rng1),
                  gt_match_get_seqid1(match),
                  rng1.start,
                  gt_range_length(&rng2),
                  gt_match_get_seqid2(match),
                  rng2.start);
  gt_file_xfputs(buffer, matchf);
}

static int process_feature(GtLTRClusterStream *lcs,
                           const char *feature,
                           GtError *err)
{
  GtArray *matches;
  GtMatchIterator *mi = NULL;
  GtMatch *match = NULL;
  GtMatchIteratorStatus status;
  char *filename,
       makeblastdb_call[BUFSIZ],
       *env,
       tmp[BUFSIZ];
  const char **argv;
  unsigned long i;
  int had_err = 0;

  if (lcs->current_state != NULL) {
    gt_free(*lcs->current_state);
    (void) snprintf(tmp, BUFSIZ, "Clustering feature: %s", feature);
    *lcs->current_state = gt_cstr_dup(tmp);
  }
  matches = gt_array_new(sizeof(GtMatch*));
  filename = (char*) gt_hashmap_get(lcs->feat_to_file, (void*) feature);
  if (!lcs->from_file) {
    argv = gt_malloc((size_t) 5 * sizeof (const char*));
    argv[0] = gt_cstr_dup("./gt suffixerator");
    argv[1] = gt_cstr_dup("-db");
    argv[2] = gt_cstr_dup(filename);
    argv[3] = gt_cstr_dup("-indexname");
    argv[4] = gt_cstr_dup(filename);
    had_err = gt_parseargsandcallsuffixerator(true, 5, argv, err);
    for (i = 0UL; i < 5UL; i++) {
      gt_free((void*) argv[i]);
    }
    gt_free((void*) argv);
  }
  if (!had_err) {
    if (!lcs->from_file) {
      /* TODO: fix problems with spaces in filename */
      env = getenv("GT_BLAST_PATH");
      if (env != NULL) {
        (void) snprintf(makeblastdb_call, BUFSIZ,
                        "%s/makeblastdb -in %s -dbtype nucl -logfile "
                        "makeblastdb.txt",
                        env, filename);
      } else {
        (void) snprintf(makeblastdb_call, BUFSIZ,
                        "makeblastdb -in %s -dbtype nucl -logfile "
                        "makeblastdb.txt",
                        filename);
      }
      had_err = system(makeblastdb_call);
    }
    if (!had_err) {
      (void) snprintf(tmp, BUFSIZ, "%s.match", filename);
      if (lcs->from_file) {
        mi = gt_match_iterator_open_new(tmp, err);
      } else
        mi = gt_match_iterator_blastn_process_new(filename,
                                                  filename,
                                                  lcs->evalue,
                                                  lcs->dust,
                                                  lcs->word_size,
                                                  lcs->gapopen,
                                                  lcs->gapextend,
                                                  lcs->penalty,
                                                  lcs->reward,
                                                  lcs->identity,
                                                  lcs->num_threads,
                                                  lcs->xdrop,
                                                  lcs->moreblast,
                                                  err);
      if (mi != NULL) {
        GtFile *matchf = NULL;
        if (!lcs->from_file)
          matchf = gt_file_new(tmp, "w", err);
        while ((status = gt_match_iterator_next(mi, &match, err))
               != GT_MATCHER_STATUS_END) {
          if (status == GT_MATCHER_STATUS_OK) {
            gt_array_add(matches, match);
            if (!lcs->from_file)
              write_match_to_file(match, matchf);
          }
        }
        gt_file_delete(matchf);
      } else
        had_err = -1;
    } else
      gt_error_set(err, "Could not run makeblastdb.");
  } else
    gt_error_set(err, "Could not run suffixerator.");
  if (!had_err) {
    GtClusteredSet *cs;
    GtEncseqLoader *encl;
    GtEncseq *encs;
    GtHashmap *seqdesc2seqnum;
    GtMatch *tmp_match;
    const char *description;
    char *output;
    unsigned long desclen,
                  num_of_seq;

    encl = gt_encseq_loader_new();
    encs = gt_encseq_loader_load(encl, filename, err);
    seqdesc2seqnum = gt_hashmap_new(GT_HASH_STRING, free_hash, NULL);
    num_of_seq = gt_encseq_num_of_sequences(encs);
    for (i = 0; i < num_of_seq; i++) {
      description = gt_encseq_description(encs, &desclen, i);
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
                            (unsigned) lcs->plarge, encs, err) != 0) {
        had_err = -1;
      }
      if (!had_err) {
        (void) cluster_annotate_nodes(cs, encs, feature, lcs->nodes, err);
      }
    } else
      had_err = -1;

    for (i = 0; i < gt_array_size(matches); i++) {
      tmp_match = *(GtMatch**) gt_array_get(matches, i);
      gt_match_delete(tmp_match);
    }
    gt_array_delete(matches);
    matches = NULL;
    gt_encseq_loader_delete(encl);
    gt_encseq_delete(encs);
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
    if (!had_err) {
      for (i = 0; i < gt_str_array_size(lcs->feat_to_file_keys); i++) {
        had_err = process_feature(lcs,
                                  gt_str_array_get(lcs->feat_to_file_keys, i),
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
  gt_hashmap_delete(lcs->feat_to_file);
  gt_str_array_delete(lcs->feat_to_file_keys);
  for (i = 0; i < gt_array_size(lcs->nodes); i++)
    gt_genome_node_delete(*(GtGenomeNode**) gt_array_get(lcs->nodes, i));
  gt_array_delete(lcs->nodes);
  gt_node_stream_delete(lcs->in_stream);
}

const GtNodeStreamClass* gt_ltr_cluster_stream_class(void)
{
  static const GtNodeStreamClass *nsc;
  if (!nsc)
    nsc = gt_node_stream_class_new(sizeof(GtLTRClusterStream),
                                   gt_ltr_cluster_stream_free,
                                   gt_ltr_cluster_stream_next);
  return nsc;
}

GtNodeStream* gt_ltr_cluster_stream_new(GtNodeStream *in_stream,
                                        GtEncseq *encseq,
                                        GtStr *file_prefix,
                                        unsigned long plarge,
                                        unsigned long psmall,
                                        double evalue,
                                        bool dust,
                                        int word_size,
                                        int gapopen,
                                        int gapextend,
                                        int penalty,
                                        int reward,
                                        int num_threads,
                                        double xdrop,
                                        double identity,
                                        const char *moreblast,
                                        bool from_file,
                                        char **current_state,
                                        GtError *err)
{
  GtNodeStream *ns;
  GtLTRClusterStream *lcs;
  ns = gt_node_stream_create(gt_ltr_cluster_stream_class(), true);
  lcs = gt_ltr_cluster_stream_cast(ns);
  lcs->in_stream = gt_node_stream_ref(in_stream);
  lcs->feat_to_file = gt_hashmap_new(GT_HASH_STRING, free_hash, free_hash);
  lcs->feat_to_file_keys = gt_str_array_new();
  lcs->nodes = gt_array_new(sizeof(GtGenomeNode*));
  lcs->lcv =
      (GtLTRClusterVisitor*) gt_ltr_cluster_visitor_new(encseq,
                                                        file_prefix,
                                                        lcs->feat_to_file,
                                                        lcs->feat_to_file_keys,
                                                        err);
  lcs->first_next = true;
  lcs->from_file = from_file;
  lcs->next_index = 0;
  lcs->plarge = plarge;
  lcs->psmall = psmall;
  lcs->evalue = evalue;
  lcs->dust = dust;
  lcs->gapopen = gapopen;
  lcs->gapextend = gapextend;
  lcs->penalty = penalty;
  lcs->reward = reward;
  lcs->num_threads = num_threads;
  lcs->word_size = word_size;
  lcs->xdrop = xdrop;
  lcs->identity = identity;
  lcs->moreblast = moreblast;
  lcs->current_state = current_state;
  return ns;
}
