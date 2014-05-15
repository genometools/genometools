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

#include <string.h>
#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "core/cstr_api.h"
#include "core/hashmap.h"
#include "core/ma_api.h"
#include "core/minmax.h"
#include "core/queue.h"
#include "core/str_api.h"
#include "core/symbol_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "extended/extract_feature_sequence.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/node_visitor_api.h"
#include "extended/region_mapping.h"
#include "extended/testcode_filter_visitor.h"

struct GtTestcodeFilterVisitor {
  const GtNodeVisitor parent_instance;
  GtRegionMapping *rm;
  GtStr *seq;
  unsigned int windowsize;
};

const GtNodeVisitorClass* gt_testcode_filter_visitor_class();

#define testcode_filter_visitor_cast(GV)\
        gt_node_visitor_cast(gt_testcode_filter_visitor_class(), GV)

/* From Fickett J.W., Nucleic Acids Res. 10(17): 5303-5318 (1982) */
static const double tcode_pos_param[4][10] =
          {{.22,.20,.34,.45,.68,.58,.93,.84,.68,.94}, /* A */
          {.23,.30,.33,.51,.48,.66,.81,.70,.70,.80},  /* C */
          {.08,.08,.16,.27,.48,.53,.64,.74,.88,.90},  /* G */
          {.09,.09,.20,.54,.44,.69,.68,.91,.97,.97}}; /* T */
static const double tcode_con_param[4][10] =
          {{.21,.81,.65,.67,.49,.62,.55,.44,.49,.28}, /* A */
          {.31,.39,.44,.43,.59,.59,.64,.51,.64,.82},  /* C */
          {.29,.33,.41,.41,.73,.64,.64,.47,.54,.40},  /* G */
          {.58,.51,.69,.56,.75,.55,.40,.39,.24,.28}}; /* T */

static inline double tcode_get_pos(double val, int pos) {
  gt_assert(pos >= 0 && pos <4);
  if (val >= 1.9) return tcode_pos_param[pos][9];
  if (val >= 1.8) return tcode_pos_param[pos][8];
  if (val >= 1.7) return tcode_pos_param[pos][7];
  if (val >= 1.6) return tcode_pos_param[pos][6];
  if (val >= 1.5) return tcode_pos_param[pos][5];
  if (val >= 1.4) return tcode_pos_param[pos][4];
  if (val >= 1.3) return tcode_pos_param[pos][3];
  if (val >= 1.2) return tcode_pos_param[pos][2];
  if (val >= 1.1) return tcode_pos_param[pos][1];
  if (val >= 0.0) return tcode_pos_param[pos][0];
  gt_assert(false); /* must be >= 0 */
}

static inline double tcode_get_con(double val, int pos) {
  gt_assert(pos >= 0 && pos <4);
  if (val >= 0.33) return tcode_con_param[pos][9];
  if (val >= 0.31) return tcode_con_param[pos][8];
  if (val >= 0.29) return tcode_con_param[pos][7];
  if (val >= 0.27) return tcode_con_param[pos][6];
  if (val >= 0.25) return tcode_con_param[pos][5];
  if (val >= 0.23) return tcode_con_param[pos][4];
  if (val >= 0.21) return tcode_con_param[pos][3];
  if (val >= 0.19) return tcode_con_param[pos][2];
  if (val >= 0.17) return tcode_con_param[pos][1];
  if (val >= 0.0) return tcode_con_param[pos][0];
  gt_assert(false); /* must be >= 0 */
}

static inline double calc_mean_tcode_window(const char *seq, int wsize) {
  double val;
  GtUword a[3] = {0,0,0},
          c[3] = {0,0,0},
          g[3] = {0,0,0},
          t[3] = {0,0,0};
  GtUword acount = 0, ccount = 0, gcount = 0, tcount = 0;
  double apos = 0.0, cpos = 0.0, gpos = 0.0, tpos = 0.0,
         afrac = 0.0, cfrac = 0.0, gfrac = 0.0, tfrac = 0.0;
  GtUword i = 0;

  for (i = 0; i < wsize; i++) {
    int codpos = i % 3;
    switch (seq[i]) {
      case 'a':
      case 'A':
        a[codpos]++; acount++;
        break;
      case 'c':
      case 'C':
        c[codpos]++; ccount++;
        break;
      case 'g':
      case 'G':
        g[codpos]++; gcount++;
        break;
      case 't':
      case 'T':
        t[codpos]++; tcount++;
        break;
      default:
        break;
    }
  }
  apos = ((double) MAX3(a[0],a[1],a[2]))/((double) MIN3(a[0],a[1],a[2])+1.0);
  cpos = ((double) MAX3(c[0],c[1],c[2]))/((double) MIN3(c[0],c[1],c[2])+1.0);
  gpos = ((double) MAX3(g[0],g[1],g[2]))/((double) MIN3(g[0],g[1],g[2])+1.0);
  tpos = ((double) MAX3(t[0],t[1],t[2]))/((double) MIN3(t[0],t[1],t[2])+1.0);
  afrac = ((double) acount)/((double) wsize);
  cfrac = ((double) ccount)/((double) wsize);
  gfrac = ((double) gcount)/((double) wsize);
  tfrac = ((double) tcount)/((double) wsize);
  val =   tcode_get_con(afrac, 0) * .11 + tcode_get_con(cfrac, 1) * .12
        + tcode_get_con(gfrac, 2) * .15 + tcode_get_con(tfrac, 3) * .14
        + tcode_get_pos(apos, 0)  * .26 + tcode_get_pos(cpos, 1)  * .18
        + tcode_get_pos(gpos, 2)  * .31 + tcode_get_pos(tpos, 3)  * .33;
  return val;
}

static double calc_mean_tcode(GtStr *seq, int wsize) {
  GtUword i = 0;
  double acc = 0.0;
  gt_assert(seq);
  for (i = 0; i < (gt_str_length(seq)-wsize)/3; i++) {
    acc += calc_mean_tcode_window(gt_str_get(seq)+(i*3), wsize);
  }
  return acc/(double)(i);
}

static int testcode_filter_visitor_feature_node(GtNodeVisitor *nv,
                                                GtFeatureNode *fn,
                                                GtError *err)
{
  GtFeatureNodeIterator *fni;
  GtFeatureNode *node;
  GtTestcodeFilterVisitor *tfv;
  tfv = testcode_filter_visitor_cast(nv);
  int had_err = 0;

  fni = gt_feature_node_iterator_new(fn);

  while (!had_err && (node = gt_feature_node_iterator_next(fni))) {
    double tcode_mean;
    gt_str_reset(tfv->seq);
    if (gt_extract_feature_sequence(tfv->seq, (GtGenomeNode*) node,
                                    gt_symbol("CDS"),
                                    false, NULL, NULL, tfv->rm, err)) {
      had_err = -1;
    }
    if (gt_str_length(tfv->seq) <= tfv->windowsize)
        continue;
    tcode_mean = calc_mean_tcode(tfv->seq, tfv->windowsize);
    gt_str_reset(tfv->seq);
    gt_str_append_double(tfv->seq, tcode_mean, 3);
    gt_feature_node_add_attribute(node, "tcode_mean", gt_str_get(tfv->seq));
  }
  gt_feature_node_iterator_delete(fni);

  return had_err;
}

void testcode_filter_visitor_free(GtNodeVisitor *nv)
{
  GtTestcodeFilterVisitor *tfv;
  if (!nv) return;
  tfv = testcode_filter_visitor_cast(nv);
  gt_region_mapping_delete(tfv->rm);
  gt_str_delete(tfv->seq);
}

const GtNodeVisitorClass* gt_testcode_filter_visitor_class()
{
  static GtNodeVisitorClass *nvc = NULL;
  gt_class_alloc_lock_enter();
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtTestcodeFilterVisitor),
                                    testcode_filter_visitor_free,
                                    NULL,
                                    testcode_filter_visitor_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}

GtNodeVisitor* gt_testcode_filter_visitor_new(GtRegionMapping *rm,
                                              unsigned int windowsize)
{
  GtTestcodeFilterVisitor *tfv;
  GtNodeVisitor *nv;
  nv = gt_node_visitor_create(gt_testcode_filter_visitor_class());
  tfv = testcode_filter_visitor_cast(nv);
  tfv->rm = gt_region_mapping_ref(rm);
  tfv->seq = gt_str_new();
  tfv->windowsize = windowsize;
  return nv;
}
