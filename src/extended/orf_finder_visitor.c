/*
  Copyright (c) 2010      Sascha Kastens <sascha.kastens@studium.uni-hamburg.de>
  Copyright (c) 2011-2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010-2012 Center for Bioinformatics, University of Hamburg

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
#include "core/codon_iterator_encseq_api.h"
#include "core/codon_iterator_simple_api.h"
#include "core/strand_api.h"
#include "core/trans_table.h"
#include "core/encseq_api.h"
#include "extended/node_visitor_api.h"
#include "extended/region_node.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/reverse_api.h"
#include "extended/orf_finder_visitor.h"
#include "extended/orf_iterator_api.h"

struct GtORFFinderVisitor {
  const GtNodeVisitor parent_instance;
  GtEncseq *encseq;
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
                        unsigned long offset, unsigned int min,
                        unsigned int max, GT_UNUSED GtError *err)
{

  gt_assert(gf);
  unsigned long tmp;

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

static int run_orffinder(GtEncseq *encseq,
                         GtFeatureNode *gf,
                         unsigned long start,
                         unsigned long end,
                         unsigned long seqid,
                         unsigned int min,
                         unsigned int max,
                         bool all,
                         GtError *err)
{
  int had_err = 0, i;
  unsigned long seqstartpos, startpos, endpos, length, sum, seqlength,
                rcstartpos;
  GtCodonIterator* ci = NULL;
  GtTranslator* translator = NULL;
  GtORFIterator* orfi = NULL;
  GtORFIteratorStatus state;
  GtRange orf_rng, tmp_orf_rng[3];
  unsigned int orf_frame;

  /* forward strand */
  seqstartpos = gt_encseq_seqstartpos(encseq, seqid);
  startpos = seqstartpos + start;
  endpos = seqstartpos + end;
  length = endpos - startpos + 1;
  ci = gt_codon_iterator_encseq_new_with_readmode(encseq, startpos, length,
                                                  GT_READMODE_FORWARD, err);
  gt_assert(ci);
  translator = gt_translator_new(ci);
  gt_assert(translator);

  orfi = gt_orf_iterator_new(ci, translator);
  gt_assert(orfi);

  for (i = 0; i < 3; i++) {
    tmp_orf_rng[i].start = GT_UNDEF_ULONG;
    tmp_orf_rng[i].end = GT_UNDEF_ULONG;
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

  if (!had_err) {
    if (!all) {
      for (i = 0; i < 3; i++) {
        if (tmp_orf_rng[i].start != GT_UNDEF_ULONG) {
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
      tmp_orf_rng[i].start = GT_UNDEF_ULONG;
      tmp_orf_rng[i].end = GT_UNDEF_ULONG;
    }

    /* reverse strand */
    if (!had_err) {
      seqlength = gt_encseq_seqlength(encseq, seqid);
      seqstartpos = (seqid == gt_encseq_num_of_sequences(encseq) -1 )
                      ? gt_encseq_total_length(encseq) - 1
                      : gt_encseq_seqstartpos(encseq, seqid+1);
      rcstartpos = gt_encseq_total_length(encseq) - 2 -
                      seqstartpos + (seqlength - end);
      ci = gt_codon_iterator_encseq_new_with_readmode(encseq, rcstartpos,
                                                      length,
                                                      GT_READMODE_REVCOMPL,
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
            if (tmp_orf_rng[i].start != GT_UNDEF_ULONG) {
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
  unsigned long seqnum = GT_UNDEF_ULONG;
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
      GtStr *val;
      val = gt_genome_node_get_seqid((GtGenomeNode*) curnode);
      if (val == NULL) {
        gt_error_set(err, "missing attribute \"seq_number\"");
        had_err = -1;
      }
      if (sscanf(gt_str_get(val), "seq%lu", &seqnum) != 1) {
        had_err = -1;
      }
      gt_assert(seqnum != GT_UNDEF_ULONG);
      if (!had_err) {
        rng = gt_genome_node_get_range((GtGenomeNode*) curnode);
        had_err = run_orffinder(lv->encseq, curnode, rng.start - 1, rng.end - 1,
                                seqnum, lv->min, lv->max, lv->all, err);
        if (gt_hashmap_get(lv->types,
                           (void*) "all") == (void*) 1) {
          break;
        }
        else if (gt_feature_node_has_children(curnode)) {
          GtFeatureNode *tmpnode = NULL;
          GtFeatureNodeIterator *tmpgfi = gt_feature_node_iterator_new(curnode);
          tmpnode = gt_feature_node_iterator_next(tmpgfi);
          while ((tmpnode = gt_feature_node_iterator_next(tmpgfi))) {
            gft = gt_feature_node_get_type(tmpnode);
            if (strcmp(gft, (const char*) GT_ORF_TYPE) == 0) {
              continue;
            }
            curnode = gt_feature_node_iterator_next(gfi);
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

GtNodeVisitor* gt_orf_finder_visitor_new(GtEncseq *encseq,
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
  lv->encseq = encseq;
  lv->types=types;
  lv->min = min;
  lv->max = max;
  lv->all = all;
  return gv;
}
