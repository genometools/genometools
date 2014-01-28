/*
  Copyright (c) 2011-2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010      Sascha Kastens <mail@skastens.de>
  Copyright (c) 2010-2013 Center for Bioinformatics, University of Hamburg

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
#include "core/class_alloc_lock.h"
#include "core/cstr_api.h"
#include "core/hashmap.h"
#include "core/ma.h"
#include "core/range.h"
#include "core/unused_api.h"
#include "core/dlist.h"
#include "core/fptr_api.h"
#include "core/hashmap.h"
#include "core/parseutils.h"
#include "core/undef_api.h"
#include "core/translator.h"
#include "core/codon_iterator_api.h"
#include "core/codon_iterator_simple_api.h"
#include "core/strand.h"
#include "core/trans_table.h"
#include "core/encseq_api.h"
#include "extended/extract_feature_sequence.h"
#include "extended/node_visitor_api.h"
#include "extended/region_node.h"
#include "extended/region_mapping.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/reverse_api.h"
#include "extended/orf_finder_visitor.h"
#include "extended/orf_iterator_api.h"

struct GtORFFinderVisitor {
  const GtNodeVisitor parent_instance;
  GtRegionMapping *rmap;
  GtHashmap *types;
  unsigned int min;
  unsigned int max;
  bool all;
};

#define GT_ORF_TYPE         "reading_frame"
#define GT_ORF_FINDER_TAG   "GenomeTools"

#define gt_orf_finder_visitor_cast(GV)\
        gt_node_visitor_cast(gt_orf_finder_visitor_class(), GV)

#define gt_orffinder_visitor_try_cast(GV)\
        gt_node_visitor_try_cast(gt_orf_finder_visitor_class(), GV)

static void orf_attach_results_to_gff3(GtFeatureNode *gf,
                                       GtRange orf_rng, unsigned int orf_frame,
                                       GtStrand strand, GT_UNUSED GtError *err)
{
  GtGenomeNode *child;
  GtStr *tag;
  tag = gt_str_new_cstr(GT_ORF_FINDER_TAG);
  if (gt_feature_node_get_strand(gf) == GT_STRAND_REVERSE)
    strand = gt_strand_invert(strand);

  orf_rng.start++; orf_rng.end++;

  GtFeatureNodeIterator *gfi;
  GtFeatureNode *curnode = NULL, *parent_node = NULL;
  GtRange gfi_range;
  char frame_buf[3];
  sprintf(frame_buf, "%d", orf_frame);

  gfi = gt_feature_node_iterator_new(gf);

  while ((curnode = gt_feature_node_iterator_next(gfi))) {
    if (strcmp(gt_feature_node_get_type(curnode),
                                              (const char*) GT_ORF_TYPE) != 0) {
      gfi_range = gt_genome_node_get_range((GtGenomeNode*) curnode);
      if (gt_range_contains(&gfi_range, &orf_rng)) {
        parent_node = curnode;
      }
    }
  }
  if (parent_node) {
    child = gt_feature_node_new(gt_genome_node_get_seqid((GtGenomeNode*) gf),
                                GT_ORF_TYPE,
                                orf_rng.start,
                                orf_rng.end,
                                strand);
    gt_feature_node_set_source((GtFeatureNode*) child, tag);
    gt_feature_node_set_attribute((GtFeatureNode*) child, "frame", frame_buf);
    gt_feature_node_add_child(parent_node,(GtFeatureNode*) child);
  }
  gt_str_delete(tag);
  gt_feature_node_iterator_delete(gfi);
}

static void process_orf(GtRange orf_rng, unsigned int orf_frame,
                        GtStrand strand, GtFeatureNode *gf,
                        GtUword offset, unsigned int min,
                        unsigned int max, GT_UNUSED GtError *err)
{

  gt_assert(gf);
  GtUword tmp;

  if ((gt_range_length(&orf_rng) >= min) &&
      (gt_range_length(&orf_rng) <= max)) {
    switch (strand) {
      case GT_STRAND_FORWARD:
        orf_rng.start = orf_rng.start + offset;
        orf_rng.end = orf_rng.end + offset;
        break;
      case GT_STRAND_REVERSE:
        tmp = orf_rng.start;
        orf_rng.start = offset - orf_rng.end;
        orf_rng.end = offset - tmp;
        break;
      default:
        exit(GT_EXIT_PROGRAMMING_ERROR);
        break;
    }
    orf_attach_results_to_gff3(gf, orf_rng, orf_frame, strand, err);
  }
}

static int run_orffinder(GtRegionMapping *rmap,
                         GtFeatureNode *gf,
                         GtUword start,
                         GT_UNUSED GtUword end,
                         unsigned int min,
                         unsigned int max,
                         bool all,
                         GtError *err)
{
  int had_err = 0, i;
  GtUword sum;
  GtCodonIterator* ci = NULL;
  GtTranslator* translator = NULL;
  GtORFIterator* orfi = NULL;
  GtORFIteratorStatus state;
  GtRange orf_rng, tmp_orf_rng[3];
  GtStr *seq = NULL;
  unsigned int orf_frame;

  /* forward strand */
  seq = gt_str_new();
  had_err = gt_extract_feature_sequence(seq,
                                        (GtGenomeNode*) gf,
                                        gt_feature_node_get_type(gf),
                                        false, NULL, NULL, rmap, err);
  if (!had_err) {
    ci = gt_codon_iterator_simple_new(gt_str_get(seq), gt_str_length(seq), err);
    gt_assert(ci);
    translator = gt_translator_new(ci);
    gt_assert(translator);

    orfi = gt_orf_iterator_new(ci, translator);
    gt_assert(orfi);

    for (i = 0; i < 3; i++) {
      tmp_orf_rng[i].start = GT_UNDEF_UWORD;
      tmp_orf_rng[i].end = GT_UNDEF_UWORD;
    }

    while ((state = gt_orf_iterator_next(orfi, &orf_rng, &orf_frame,
                                                err)) == GT_ORF_ITERATOR_OK) {
        if (all) {
          process_orf(orf_rng, orf_frame, GT_STRAND_FORWARD, gf,
                      start, min, max, err);
        } else {
          if (gt_range_length(&orf_rng) >
              gt_range_length(&tmp_orf_rng[orf_frame])) {
            tmp_orf_rng[orf_frame].start = orf_rng.start;
            tmp_orf_rng[orf_frame].end = orf_rng.end;
          }
        }
    }
    if (state == GT_ORF_ITERATOR_ERROR)
      had_err = -1;
  }

  if (!had_err) {
    if (!all) {
      for (i = 0; i < 3; i++) {
        if (tmp_orf_rng[i].start != GT_UNDEF_UWORD) {
          process_orf(tmp_orf_rng[i], (unsigned int) i, GT_STRAND_FORWARD, gf,
                      start, min, max, err);
        }
      }
    }
    gt_codon_iterator_delete(ci);
    gt_translator_delete(translator);
    gt_orf_iterator_delete(orfi);
    orfi = NULL;
    ci = NULL;
    translator = NULL;

    for (i = 0; i < 3; i++) {
      tmp_orf_rng[i].start = GT_UNDEF_UWORD;
      tmp_orf_rng[i].end = GT_UNDEF_UWORD;
    }

    /* reverse strand */
    if (!had_err) {
      GT_UNUSED int rval = 0;
      GtUword length = gt_str_length(seq);
      char *strp = (char*) gt_str_get_mem(seq);
      rval = gt_reverse_complement(strp, gt_str_length(seq), err);
      gt_assert(!rval); /* XXX */
      ci = gt_codon_iterator_simple_new(gt_str_get(seq), gt_str_length(seq),
                                        err);
      gt_assert(ci);
      translator = gt_translator_new(ci);
      gt_assert(translator);
      orfi = gt_orf_iterator_new(ci, translator);
      gt_assert(orfi);

      sum = start + length - 1;

      while ((state = gt_orf_iterator_next(orfi, &orf_rng, &orf_frame,
                                                  err)) == GT_ORF_ITERATOR_OK) {
          if (all) {
            process_orf(orf_rng, orf_frame, GT_STRAND_REVERSE, gf,
                        sum, min, max, err);
          } else {
            if (gt_range_length(&orf_rng) >
                gt_range_length(&tmp_orf_rng[orf_frame])) {
              tmp_orf_rng[orf_frame].start = orf_rng.start;
              tmp_orf_rng[orf_frame].end = orf_rng.end;
            }
          }
      }
      if (state == GT_ORF_ITERATOR_ERROR)
        had_err = -1;
      if (!had_err) {
        if (!all) {
          for (i = 0; i < 3; i++) {
            if (tmp_orf_rng[i].start != GT_UNDEF_UWORD) {
              process_orf(tmp_orf_rng[i], (unsigned int) i, GT_STRAND_REVERSE,
                          gf, sum, min, max, err);
            }
          }
        }
      }
    }
    gt_codon_iterator_delete(ci);
    gt_translator_delete(translator);
    gt_orf_iterator_delete(orfi);
  }
  gt_str_delete(seq);
  return had_err;
}

static int gt_orf_finder_visitor_feature_node(GtNodeVisitor *gv,
                                              GtFeatureNode *gf,
                                              GtError *err)
{
  GtORFFinderVisitor *lv;
  const char *gft = NULL;
  GtFeatureNodeIterator *gfi;
  GtFeatureNode *curnode = NULL;
  int had_err = 0;
  GtRange rng;

  lv = gt_orf_finder_visitor_cast(gv);
  gt_assert(lv);
  gt_error_check(err);

  gfi = gt_feature_node_iterator_new(gf);

  while (!had_err && (curnode = gt_feature_node_iterator_next(gfi))) {
    gft = gt_feature_node_get_type(curnode);

    if (gt_hashmap_get(lv->types, (void*) gft) != NULL ||
                       gt_hashmap_get(lv->types,
                                      (void*) "all") == (void*) 1) {
      if (!had_err) {
        rng = gt_genome_node_get_range((GtGenomeNode*) curnode);
        had_err = run_orffinder(lv->rmap, curnode, rng.start - 1, rng.end - 1,
                                lv->min, lv->max, lv->all, err);
      }
      if (!had_err) {
        if (gt_hashmap_get(lv->types,
                           (void*) "all") == (void*) 1) {
          break;
        }
        else if (gt_feature_node_has_children(curnode)) {
          GtFeatureNode *tmpnode = NULL;
          GtFeatureNodeIterator *tmpgfi = gt_feature_node_iterator_new(curnode);
          (void) gt_feature_node_iterator_next(tmpgfi);
          while ((tmpnode = gt_feature_node_iterator_next(tmpgfi))) {
            gft = gt_feature_node_get_type(tmpnode);
            if (strcmp(gft, (const char*) GT_ORF_TYPE) == 0) {
              continue;
            }
            /* curnode = gt_feature_node_iterator_next(gfi); */
          }
          gt_feature_node_iterator_delete(tmpgfi);
        }
      }
    }
  }

  gt_feature_node_iterator_delete(gfi);

  return had_err;
}

const GtNodeVisitorClass* gt_orf_finder_visitor_class(void)
{
  static const GtNodeVisitorClass *gvc = NULL;
  gt_class_alloc_lock_enter();
  if (!gvc) {
    gvc = gt_node_visitor_class_new(sizeof (GtORFFinderVisitor),
                                    NULL,
                                    NULL,
                                    gt_orf_finder_visitor_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  gt_class_alloc_lock_leave();
  return gvc;
}

GtNodeVisitor* gt_orf_finder_visitor_new(GtRegionMapping *rmap,
                                         GtHashmap *types,
                                         unsigned int min,
                                         unsigned int max,
                                         bool all,
                                         GT_UNUSED GtError *err)
{
  GtNodeVisitor *gv;
  GtORFFinderVisitor *lv;
  gv = gt_node_visitor_create(gt_orf_finder_visitor_class());
  lv = gt_orf_finder_visitor_cast(gv);
  gt_assert(lv);
  lv->rmap = rmap;
  lv->types=types;
  lv->min = min;
  lv->max = max;
  lv->all = all;
  return gv;
}
