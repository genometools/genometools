/*
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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
#include "core/array_api.h"
#include "core/class_alloc_lock.h"
#include "core/cstr_api.h"
#include "core/encseq_api.h"
#include "core/ma.h"
#include "core/parseutils.h"
#include "core/range.h"
#include "core/str_api.h"
#include "core/unused_api.h"
#include "extended/node_visitor_api.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/feature_type.h"
#include "ltr/ltrharvest_tabout_visitor.h"

struct GtLTRharvestTaboutVisitor {
  const GtNodeVisitor parent_instance;
  bool longoutput;
  const GtEncseq *encseq;
};

#define gt_ltrharvest_tabout_visitor_cast(GV)\
        gt_node_visitor_cast(gt_ltrharvest_tabout_visitor_class(), GV)

static inline void
gt_ltrharvest_tabout_visitor_seq4feat(GtLTRharvestTaboutVisitor *nv,
                                      GtFeatureNode *fn, GtStr *out,
                                      unsigned long seqnum)
{
  char *buf;
  GtRange rng;
  unsigned long startpos;
  gt_assert(nv && fn && out);
  rng = gt_genome_node_get_range((GtGenomeNode*) fn);
  buf = gt_calloc((size_t) (gt_range_length(&rng)+1), sizeof (char));
  startpos = gt_encseq_seqstartpos(nv->encseq, seqnum);
  gt_encseq_extract_decoded(nv->encseq, buf,
                            startpos + rng.start - 1,
                            startpos + rng.end - 1);
  gt_str_append_cstr(out, buf);
  gt_free(buf);
}

static int gt_ltrharvest_tabout_visitor_feature_node(GtNodeVisitor *nv,
                                                     GtFeatureNode *fn,
                                                     GT_UNUSED GtError *err)
{
  GtLTRharvestTaboutVisitor *lv;
  GtFeatureNodeIterator *fni;
  GtFeatureNode *curnode = NULL,
                *leftltr = NULL,
                *rightltr = NULL,
                *ltr_retrotrans = NULL,
                *lefttsd = NULL,
                *righttsd = NULL,
                *leftmotif_5 = NULL,
                *leftmotif_3 = NULL,
                *rightmotif_5 = NULL,
                *rightmotif_3 = NULL;
  int had_err = 0,
      added_motifs = 0,
      added_tsd = 0,
      added_ltr = 0;
  GtStr *line;
  GtRange rng;
  bool no_element = false;
  unsigned long seqnum = 0;
  const char *fnt;
  char buf[BUFSIZ];
  lv = gt_ltrharvest_tabout_visitor_cast(nv);
  gt_assert(lv);
  gt_error_check(err);

  line = gt_str_new();
  /* traverse annotation subgraph and collect info */
  fni = gt_feature_node_iterator_new(fn);
  while (!had_err && (curnode = gt_feature_node_iterator_next(fni))) {
    fnt = gt_feature_node_get_type(curnode);

    if (strcmp(fnt, gt_ft_LTR_retrotransposon) == 0)
    {
      const char *val;
      ltr_retrotrans = curnode;
      val = gt_feature_node_get_attribute(curnode, "seq_number");
      if (val == NULL) {
        gt_error_set(err, "missing attribute \"seq_number\"");
        had_err = -1;
      }
      if (!had_err) {
        (void) gt_parse_ulong(&seqnum, val);
      }
    }
    if (strcmp(fnt, gt_ft_long_terminal_repeat) == 0)
    {
      switch (added_ltr) {
        case 0:
          leftltr = curnode;
          break;
        case 1:
          rightltr = curnode;
          break;
        default:
          gt_error_set(err,
                       "invalid number of LTR annotations: more than 2");
          had_err = -1;
      }
      added_ltr++;
    }
    if (strcmp(fnt, gt_ft_inverted_repeat) == 0)
    {
      switch (added_motifs) {
        case 0:
          leftmotif_5 = curnode;
          break;
        case 1:
          leftmotif_3 = curnode;
          break;
        case 2:
          rightmotif_5 = curnode;
          break;
        case 3:
          rightmotif_3 = curnode;
          break;
        default:
          gt_error_set(err,
                       "invalid number of motif annotations: more than 4");
          had_err = -1;
          break;
      }
      added_motifs++;
    }
    if (strcmp(fnt, gt_ft_target_site_duplication) == 0)
    {
      switch (added_tsd) {
        case 0:
          lefttsd = curnode;
          break;
        case 1:
          righttsd = curnode;
          break;
        default:
          gt_error_set(err,
                       "invalid number of TSD annotations: more than 2");
          had_err = -1;
      }
      added_tsd++;
    }
  }
  gt_feature_node_iterator_delete(fni);

  /* check for invalid annotations */
  if (!had_err && (!leftltr || !rightltr)) {
    gt_error_set(err, "missing LTR annotations");
    had_err = -1;
  }
  if (lv->longoutput) {
    if (!had_err && added_motifs > 0 && added_motifs != 4) {
      gt_error_set(err, "invalid number of motif annotations, != 4");
      had_err = -1;
    }
    if (!had_err && added_tsd > 0 && added_tsd != 2) {
      gt_error_set(err, "invalid number of TSD annotations, != 2");
      had_err = -1;
    }
  }

  /* output in tabular format */
  if (!had_err) {
    if (ltr_retrotrans != NULL) {
      /* whole element */
      rng = gt_genome_node_get_range((GtGenomeNode*) ltr_retrotrans);
      (void) snprintf(buf, BUFSIZ-1, "%lu  %lu  %lu  ", rng.start,
                      rng.end, gt_range_length(&rng));
      gt_str_append_cstr(line, buf);
    } else {
      no_element = true;
    }
  }

  if (!had_err && !no_element) {
    /* left LTR */
    rng = gt_genome_node_get_range((GtGenomeNode*) leftltr);
    (void) snprintf(buf, BUFSIZ-1, "%lu  %lu  %lu  ", rng.start,
                    rng.end, gt_range_length(&rng));
    gt_str_append_cstr(line, buf);
    /* left TSD */
    if (lv->longoutput && added_tsd > 0 ) {
      gt_assert(lefttsd);
      gt_ltrharvest_tabout_visitor_seq4feat(lv, lefttsd, line, seqnum);
      gt_str_append_cstr(line, "  ");
      rng = gt_genome_node_get_range((GtGenomeNode*) lefttsd);
      gt_str_append_ulong(line, gt_range_length(&rng));
      gt_str_append_cstr(line, "  ");
    }
    /* left motif */
    if (lv->longoutput && added_motifs > 0 ) {
      gt_assert(leftmotif_5 && leftmotif_3);
      gt_ltrharvest_tabout_visitor_seq4feat(lv, leftmotif_5, line, seqnum);
      gt_str_append_cstr(line, "..");
      gt_ltrharvest_tabout_visitor_seq4feat(lv, leftmotif_3, line, seqnum);
      gt_str_append_cstr(line, "  ");
    }

    /* right LTR */
    rng = gt_genome_node_get_range((GtGenomeNode*) rightltr);
    (void) snprintf(buf, BUFSIZ-1, "%lu  %lu  %lu  ", rng.start,
                    rng.end, gt_range_length(&rng));
    gt_str_append_cstr(line, buf);
    /* right TSD */
    if (lv->longoutput && added_tsd > 0 ) {
      gt_assert(righttsd);
      gt_ltrharvest_tabout_visitor_seq4feat(lv, righttsd, line, seqnum);
      gt_str_append_cstr(line, "  ");
      rng = gt_genome_node_get_range((GtGenomeNode*) righttsd);
      gt_str_append_ulong(line, gt_range_length(&rng));
      gt_str_append_cstr(line, "  ");
    }
    /* right motif */
    if (lv->longoutput && added_motifs > 0 ) {
      gt_assert(rightmotif_5 && rightmotif_3);
      gt_ltrharvest_tabout_visitor_seq4feat(lv, rightmotif_5, line, seqnum);
      gt_str_append_cstr(line, "..");
      gt_ltrharvest_tabout_visitor_seq4feat(lv, rightmotif_3, line, seqnum);
      gt_str_append_cstr(line, "  ");
    }

    gt_str_append_cstr(line, gt_feature_node_get_attribute(ltr_retrotrans,
                                                           "ltr_similarity"));
    gt_str_append_cstr(line, "  ");
    gt_str_append_ulong(line, seqnum);
    printf("%s\n", gt_str_get(line));
  }
  gt_str_delete(line);
  return had_err;
}

const GtNodeVisitorClass* gt_ltrharvest_tabout_visitor_class(void)
{
  static const GtNodeVisitorClass *nvc = NULL;
  gt_class_alloc_lock_enter();
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtLTRharvestTaboutVisitor),
                                    NULL,
                                    NULL,
                                    gt_ltrharvest_tabout_visitor_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}

GtNodeVisitor* gt_ltrharvest_tabout_visitor_new(void)
{
  GtNodeVisitor *nv;
  GtLTRharvestTaboutVisitor *lv;
  nv = gt_node_visitor_create(gt_ltrharvest_tabout_visitor_class());
  lv = gt_ltrharvest_tabout_visitor_cast(nv);
  gt_assert(lv);
  lv->longoutput = false;
  lv->encseq = NULL;
  return nv;
}

GtNodeVisitor* gt_ltrharvest_tabout_visitor_new_longoutput(const GtEncseq
                                                                        *encseq)
{
  GtNodeVisitor *nv;
  GtLTRharvestTaboutVisitor *lv;
  gt_assert(encseq);
  nv = gt_node_visitor_create(gt_ltrharvest_tabout_visitor_class());
  lv = gt_ltrharvest_tabout_visitor_cast(nv);
  gt_assert(lv);
  lv->longoutput = true;
  lv->encseq = encseq;
  return nv;
}
