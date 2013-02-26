/*
  Copyright (c) 2012-2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2012-2013 Center for Bioinformatics, University of Hamburg

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
#include "extended/snp_annotator_stream.h"
#include "extended/snp_annotator_visitor.h"

  struct GtSNPAnnotatorStream {
  const GtNodeStream parent_instance;
  GtNodeStream *merge_stream;
  GtTransTable *tt;
  GtArray *cur_gene_set,
          *instreams;
  GtRange cur_gene_range;
  GtRegionMapping *rmap;
  GtQueue *snps,
          *outqueue;
};

#define gt_snp_annotator_stream_cast(GS)\
        gt_node_stream_cast(gt_snp_annotator_stream_class(), GS);

static int snp_annotator_stream_process_snp(void **elem, void *info,
                                            GtError *err)
{
  int had_err = 0;
  GtGenomeNode *snp = *(GtGenomeNode**) elem;
  GtNodeVisitor *sav = (GtNodeVisitor*) info;
  had_err = gt_genome_node_accept(snp, sav, err);
  return had_err;
}

static int snp_annotator_stream_process_current_gene(GtSNPAnnotatorStream *sas,
                                                     GtError *err)
{
  int had_err = 0;
  unsigned long i;
  unsigned long nof_genes = gt_array_size(sas->cur_gene_set);
  gt_error_check(err);

  if (gt_queue_size(sas->snps) > 0) {
    /* we need to process SNPs for a gene cluster*/
    gt_assert(gt_queue_size(sas->outqueue) == 0);

    for (i = 0; !had_err && i < nof_genes; i++) {
      GtNodeVisitor *sav;
      GtFeatureNode *gene;
      gene = *(GtFeatureNode**) gt_array_get(sas->cur_gene_set, i);
      sav = gt_snp_annotator_visitor_new(gene, sas->tt, sas->rmap, err);
      if (!sav)
        had_err = -1;
      if (!had_err) {
        if (i < nof_genes-1) {
          had_err = gt_queue_iterate(sas->snps,
                                     snp_annotator_stream_process_snp,
                                     sav, err);
        } else {
          while (!had_err && gt_queue_size(sas->snps) > 0) {
            GtFeatureNode *snp = (GtFeatureNode*) gt_queue_get(sas->snps);
            had_err = gt_genome_node_accept((GtGenomeNode*) snp, sav, err);
            gt_queue_add(sas->outqueue, snp);
            gt_genome_node_delete((GtGenomeNode*) snp);
          }
        }
        gt_node_visitor_delete(sav);
      }
      gt_genome_node_delete((GtGenomeNode*) gene);
    }
  } else {
    /* no SNPs for this gene cluster, delete it */
    for (i = 0; !had_err && i < nof_genes; i++) {
      gt_genome_node_delete(*(GtGenomeNode**) gt_array_get(sas->cur_gene_set,
                                                           i));
    }
  }
  gt_assert(gt_queue_size(sas->snps) == 0);
  gt_array_reset(sas->cur_gene_set);
  return had_err;
}

static int snp_annotator_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                     GtError *err)
{
  GtSNPAnnotatorStream *sas;
  int had_err = 0;
  bool complete_cluster = false;
  GtGenomeNode *mygn = NULL;
  GtFeatureNode *fn = NULL;
  const char *snv_type = gt_symbol(gt_ft_SNV),
             *snp_type = gt_symbol(gt_ft_SNP),
             *gene_type = gt_symbol(gt_ft_gene);
  gt_error_check(err);
  sas = gt_snp_annotator_stream_cast(ns);

  /* if there are still SNPs left in the buffer, output them */
  if (gt_queue_size(sas->outqueue) > 0) {
    *gn = (GtGenomeNode*) gt_queue_get(sas->outqueue);
    return had_err;
  } else complete_cluster = false;

  while (!had_err && !complete_cluster) {
    had_err = gt_node_stream_next(sas->merge_stream, &mygn, err);

    /* stop if stream is at the end */
    if (had_err || !mygn) break;

    /* process all feature nodes */
    if ((fn = gt_feature_node_try_cast(mygn))) {
      GtGenomeNode *addgn;
      const char *type = gt_feature_node_get_type(fn);
      GtRange new_rng = gt_genome_node_get_range(mygn);
      if (type == snv_type || type == snp_type) {
        /* -----> this is a SNP <----- */
        if (gt_range_overlap(&new_rng, &sas->cur_gene_range)) {
          /* it falls into the currently observed range */
          gt_queue_add(sas->snps, gt_genome_node_ref((GtGenomeNode*) fn));
        } else {
          /* SNP outside a gene, this cluster is done
             add to out queue and start serving */
          gt_assert(gt_queue_size(sas->outqueue) == 0);
          had_err = snp_annotator_stream_process_current_gene(sas, err);
          gt_queue_add(sas->outqueue, mygn);
          if (gt_queue_size(sas->outqueue) > 0) {
            *gn = (GtGenomeNode*) gt_queue_get(sas->outqueue);
            complete_cluster = true;
          }
        }
      } else if (type == gene_type) {
        /* -----> this is a gene <----- */
        if (gt_array_size(sas->cur_gene_set) == 0UL) {
          /* new overlapping gene cluster */
          addgn = gt_genome_node_ref(mygn);
          gt_array_add(sas->cur_gene_set, addgn);
          sas->cur_gene_range = gt_genome_node_get_range(mygn);
        } else {
          if (gt_range_overlap(&new_rng, &sas->cur_gene_range)) {
            /* gene overlaps with current one, add to cluster */
            addgn = gt_genome_node_ref(mygn);
            gt_array_add(sas->cur_gene_set, addgn);
            sas->cur_gene_range = gt_range_join(&sas->cur_gene_range, &new_rng);
          } else {
            /* finish current cluster and start a new one */
            had_err = snp_annotator_stream_process_current_gene(sas, err);
            if (!had_err) {
              addgn = gt_genome_node_ref(mygn);
              gt_array_add(sas->cur_gene_set, addgn);
              sas->cur_gene_range = gt_genome_node_get_range(mygn);
            }
            if (gt_queue_size(sas->outqueue) > 0) {
              *gn = (GtGenomeNode*) gt_queue_get(sas->outqueue);
              complete_cluster = true;
            }
          }
        }
        /* from now on, genes are kept in gene cluster arrays only */
        gt_genome_node_delete(mygn);
      }
    } else {
      /* meta node */
      had_err = snp_annotator_stream_process_current_gene(sas, err);
      if (!had_err) {
        gt_queue_add(sas->outqueue, mygn);
      }
      if (gt_queue_size(sas->outqueue) > 0) {
        *gn = (GtGenomeNode*) gt_queue_get(sas->outqueue);
        complete_cluster = true;
      }
    }
  }

  return had_err;
}

static void snp_annotator_stream_free(GtNodeStream *ns)
{
  unsigned long i;
  GtSNPAnnotatorStream *sas;
  if (!ns) return;
  sas = gt_snp_annotator_stream_cast(ns);
  gt_region_mapping_delete(sas->rmap);
  while (gt_queue_size(sas->snps) > 0) {
    gt_genome_node_delete((GtGenomeNode*) gt_queue_get(sas->snps));
  }
  while (gt_queue_size(sas->outqueue) > 0) {
    gt_genome_node_delete((GtGenomeNode*) gt_queue_get(sas->outqueue));
  }
  for (i = 0; i < gt_array_size(sas->instreams); i++) {
    gt_node_stream_delete(*(GtNodeStream**) gt_array_get(sas->instreams, i));
  }
  for (i = 0; i < gt_array_size(sas->cur_gene_set); i++) {
    gt_genome_node_delete(*(GtGenomeNode**) gt_array_get(sas->cur_gene_set, i));
  }
  gt_array_delete(sas->cur_gene_set);
  gt_node_stream_delete(sas->merge_stream);
  gt_array_delete(sas->instreams);
  gt_queue_delete(sas->snps);
  gt_queue_delete(sas->outqueue);
}

const GtNodeStreamClass* gt_snp_annotator_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtSNPAnnotatorStream),
                                   snp_annotator_stream_free,
                                   snp_annotator_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

GtNodeStream* gt_snp_annotator_stream_new(GtNodeStream *gvf_stream,
                                          GtNodeStream *gff_stream,
                                          GtTransTable *trans_table,
                                          GtRegionMapping *rmap)
{
  GtSNPAnnotatorStream *sas;
  GtNodeStream *ns;
  gt_assert(gvf_stream && gff_stream && rmap);
  ns = gt_node_stream_create(gt_snp_annotator_stream_class(), true);
  sas = gt_snp_annotator_stream_cast(ns);
  sas->instreams = gt_array_new(sizeof (GtNodeStream*));
  (void) gt_node_stream_ref(gvf_stream);
  gt_array_add(sas->instreams, gvf_stream);
  (void) gt_node_stream_ref(gff_stream);
  gt_array_add(sas->instreams, gff_stream);
  sas->cur_gene_set = gt_array_new(sizeof (GtFeatureNode*));
  sas->merge_stream = gt_merge_stream_new(sas->instreams);
  sas->rmap = gt_region_mapping_ref(rmap);
  sas->cur_gene_range.start =
    sas->cur_gene_range.end = GT_UNDEF_ULONG;
  sas->snps = gt_queue_new();
  sas->outqueue = gt_queue_new();
  sas->tt = trans_table;
  return ns;
}
