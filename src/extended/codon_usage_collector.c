/*
  Copyright (c) 2014 Sascha Steinbiss <ss34@sanger.ac.uk>

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

#include <math.h>
#include <string.h>
#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "core/codon_iterator_api.h"
#include "core/codon_iterator_simple_api.h"
#include "core/ma_api.h"
#include "core/str_api.h"
#include "core/symbol_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/codon_usage_collector.h"
#include "extended/extract_feature_sequence.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/node_visitor_api.h"
#include "extended/region_mapping.h"

struct GtCodonUsageCollector {
  const GtNodeVisitor parent_instance;
  GtRegionMapping *rm;
  GtStr *seq;
  GtUword codon_count[64],
          k;
  double codon_freqs[64];
};

const GtNodeVisitorClass* gt_codon_usage_collector_class();

#define codon_usage_collector_cast(GV)\
        gt_node_visitor_cast(gt_codon_usage_collector_class(), GV)

#define CUC_A_CODE  0
#define CUC_C_CODE  1
#define CUC_G_CODE  2
#define CUC_T_CODE  3

static inline int codon_usage_collector_codon2idx(char c0, char c1, char c2)
{
  unsigned int code = 0;
  switch (c0)
  {
    case 'a':
    case 'A':
        code = (CUC_A_CODE << 4);
        break;
    case 'c':
    case 'C':
        code = (CUC_C_CODE << 4);
        break;
    case 'g':
    case 'G':
        code = (CUC_G_CODE << 4);
        break;
    case 't':
    case 'T':
        code = (CUC_T_CODE << 4);;
        break;
    default:
        return GT_UNDEF_INT;
  }
  switch (c1)
  {
    case 'a':
    case 'A':
        code += (CUC_A_CODE << 2);
        break;
    case 'c':
    case 'C':
        code += (CUC_C_CODE << 2);
        break;
    case 'g':
    case 'G':
        code += (CUC_G_CODE << 2);
        break;
    case 't':
    case 'T':
        code += (CUC_T_CODE << 2);
        break;
    default:
        return GT_UNDEF_INT;
  }
  switch (c2)
  {
    case 'a':
    case 'A':
        code += CUC_A_CODE;
        break;
    case 'c':
    case 'C':
        code += CUC_C_CODE;
        break;
    case 'g':
    case 'G':
        code += CUC_G_CODE;
        break;
    case 't':
    case 'T':
        code += CUC_T_CODE;
        break;
    default:
        return GT_UNDEF_INT;
  }
  return code;
}

static int codon_usage_collector_calc_codon_usage(GtCodonUsageCollector *cuc,
                                                  GtStr *seq, GtError *err)
{
  int had_err = 0, rval = 0;
  GT_UNUSED unsigned int frame;
  char c1, c2, c3;
  GtCodonIterator *ci = gt_codon_iterator_simple_new(gt_str_get(seq),
                                                     gt_str_length(seq), err);
  if (!ci)
    return -1;

  while (!(rval = gt_codon_iterator_next(ci, &c1, &c2, &c3, &frame, err))) {
    int code = codon_usage_collector_codon2idx(c1, c2, c3);
    if (code != GT_UNDEF_INT && frame == 0) {
      gt_assert(code < 64);
      cuc->codon_count[code]++;
      cuc->k++;
    }
  }

  if (rval == GT_CODON_ITERATOR_ERROR)
    had_err = -1;

  gt_codon_iterator_delete(ci);

  return had_err;
}

static int codon_usage_collector_feature_node(GtNodeVisitor *nv,
                                              GtFeatureNode *fn,
                                              GtError *err)
{
  GtFeatureNodeIterator *fni;
  GtFeatureNode *node;
  GtCodonUsageCollector *cuc;
  cuc = codon_usage_collector_cast(nv);
  int had_err = 0;

  fni = gt_feature_node_iterator_new(fn);
  while (!had_err && (node = gt_feature_node_iterator_next(fni))) {
    gt_str_reset(cuc->seq);
    if (gt_feature_node_get_type(node) != gt_symbol("CDS"))
      continue;
    if (gt_extract_feature_sequence(cuc->seq, (GtGenomeNode*) node,
                                    gt_symbol("CDS"),
                                    false, NULL, NULL, cuc->rm, err)) {
      had_err = -1;
    }
    if (!had_err)
      had_err = codon_usage_collector_calc_codon_usage(cuc, cuc->seq, err);
  }
  gt_feature_node_iterator_delete(fni);

  return had_err;
}

void codon_usage_collector_free(GtNodeVisitor *nv)
{
  GtCodonUsageCollector *cuc;
  if (!nv) return;
  cuc = codon_usage_collector_cast(nv);
  gt_region_mapping_delete(cuc->rm);
  gt_str_delete(cuc->seq);
}

const GtNodeVisitorClass* gt_codon_usage_collector_class()
{
  static GtNodeVisitorClass *nvc = NULL;
  gt_class_alloc_lock_enter();
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtCodonUsageCollector),
                                    codon_usage_collector_free,
                                    NULL,
                                    codon_usage_collector_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}

GtNodeVisitor* gt_codon_usage_collector_new(GtRegionMapping *rm)
{
  GtCodonUsageCollector *cuc;
  GtNodeVisitor *nv;
  nv = gt_node_visitor_create(gt_codon_usage_collector_class());
  cuc = codon_usage_collector_cast(nv);
  cuc->rm = gt_region_mapping_ref(rm);
  cuc->seq = gt_str_new();
  memset(cuc->codon_count, 0, 64 * (sizeof (GtUword)));
  cuc->k = 0;
  return nv;
}

void gt_codon_usage_collector_output_usage_table(GtCodonUsageCollector *cuc,
                                                 GtFile *outfile)
{
  int i = 0;
  double acc = 0.0, fmeans;
  GT_UNUSED GtUword tmp = 0;
  gt_assert(cuc);

  if (cuc->k == 0)
    return;

  /* check whether the counts match up with total codons */
  for (i=0; i < 64; i++) {
    tmp += cuc->codon_count[i];
  }
  gt_assert(tmp == cuc->k);

  /* logify, adjust, output */
  for (i=0; i < 64; i++) {
    if (cuc->codon_count[i] == 0) {
      cuc->codon_freqs[i] = log(1.0/((double)cuc->k));
    } else {
      cuc->codon_freqs[i] =
                           log((double) cuc->codon_count[i]/((double)cuc->k));
    }
    acc += cuc->codon_freqs[i];
  }
  fmeans = acc/64;
  /* XXX use GtTransTable to allow for different stop codons */
  cuc->codon_freqs[codon_usage_collector_codon2idx('t','a','g')] = fmeans;
  cuc->codon_freqs[codon_usage_collector_codon2idx('t','a','a')] = fmeans;
  cuc->codon_freqs[codon_usage_collector_codon2idx('t','g','a')] = fmeans;

  for (i=0; i < 64; i++) {
    gt_file_xprintf(outfile, GT_WU "\t" GT_WU "\t%f\n",
                    cuc->codon_count[i],
                    cuc->k,
                    cuc->codon_freqs[i]);
  }
}
