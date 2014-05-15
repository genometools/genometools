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

#include <math.h>
#include <string.h>
#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "core/codon_iterator_api.h"
#include "core/codon_iterator_simple_api.h"
#include "core/file.h"
#include "core/log_api.h"
#include "core/ma_api.h"
#include "core/minmax.h"
#include "core/parseutils_api.h"
#include "core/str_api.h"
#include "core/strand_api.h"
#include "core/symbol_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "extended/extract_feature_sequence.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/node_visitor_api.h"
#include "extended/region_mapping.h"
#include "extended/codon_usage_visitor.h"

struct GtCodonUsageVisitor {
  const GtNodeVisitor parent_instance;
  GtRegionMapping *rm;
  GtStr *seq;
  GtUword windowsize, startpos;
  double codon_freqs[64], coding_threshold;
};

const GtNodeVisitorClass* gt_codon_usage_visitor_class();

#define codon_usage_visitor_cast(GV)\
        gt_node_visitor_cast(gt_codon_usage_visitor_class(), GV)

#define CUV_A_CODE  0
#define CUV_C_CODE  1
#define CUV_G_CODE  2
#define CUV_T_CODE  3

static inline int codon_usage_visitor_codon2idx(char c0, char c1, char c2)
{
  unsigned int code = 0;
  switch (c0)
  {
    case 'a':
    case 'A':
        code = (CUV_A_CODE << 4);
        break;
    case 'c':
    case 'C':
        code = (CUV_C_CODE << 4);
        break;
    case 'g':
    case 'G':
        code = (CUV_G_CODE << 4);
        break;
    case 't':
    case 'T':
        code = (CUV_T_CODE << 4);;
        break;
    default:
        return GT_UNDEF_INT;
  }
  switch (c1)
  {
    case 'a':
    case 'A':
        code += (CUV_A_CODE << 2);
        break;
    case 'c':
    case 'C':
        code += (CUV_C_CODE << 2);
        break;
    case 'g':
    case 'G':
        code += (CUV_G_CODE << 2);
        break;
    case 't':
    case 'T':
        code += (CUV_T_CODE << 2);
        break;
    default:
        return GT_UNDEF_INT;
  }
  switch (c2)
  {
    case 'a':
    case 'A':
        code += CUV_A_CODE;
        break;
    case 'c':
    case 'C':
        code += CUV_C_CODE;
        break;
    case 'g':
    case 'G':
        code += CUV_G_CODE;
        break;
    case 't':
    case 'T':
        code += CUV_T_CODE;
        break;
    default:
        return GT_UNDEF_INT;
  }
  return code;
}

static inline int codon_usage_visitor_calc_window(GtCodonUsageVisitor *cuv,
                                                  const char *seq,
                                                  double *rval,
                                                  GtError *err)
{
  int had_err = 0, status;
  unsigned int frame;
  char c1, c2, c3;
  double h[3],
         p[3],
         logp[3],
         mean,
         q = (1.0/3.0),
         denom;
  GtUword i = 0, nof_codons = 0;
  GtCodonIterator *ci;

  ci = gt_codon_iterator_simple_new(seq, cuv->windowsize, err);
  if (!ci)
    return -1;
  h[0] = h[1] = h[2] = 0.0;
  while (!(status = gt_codon_iterator_next(ci, &c1, &c2, &c3, &frame, err))) {
    int code = codon_usage_visitor_codon2idx(c1, c2, c3);
    if (code != GT_UNDEF_INT) {
      gt_assert(code < 64);
      h[(frame + cuv->startpos) % 3] += cuv->codon_freqs[code];
      nof_codons++;
    }
  }
  gt_codon_iterator_delete(ci);
  if (status == GT_CODON_ITERATOR_ERROR)
    had_err = -1;

  if (!had_err) {
    mean = (h[0] + h[1] + h[2])/3;
    for (i = 0; i < 3; i++) {
      h[i] -= mean;
    }
    denom = (q*exp(h[0]))+(q*exp(h[1]))+(q*exp(h[2]));
    gt_assert(denom != 0.0);
    for (i = 0; i < 3; i++) {
      p[i] = (q*exp(h[i]))/((double) denom);
      if (p[i] == 1) {
        /* XXX: hack */
        p[i] = 0.999999;
      }
      gt_assert(!isinf(p[i]));
      logp[i] = log10(p[i]/(1.0-p[i]));
      rval[i] = logp[i];
    }
  }

  return had_err;
}

static int codon_usage_visitor_calc_value(GtCodonUsageVisitor *cuv,
                                          GtStr *seq,
                                          double *result,
                                          GtError *err) {
  GtUword i = 0, codons;
  double acc[3] = {0.0,0.0,0.0};
  int had_err = 0;
  gt_assert(seq && result);

  codons = (gt_str_length(seq) - cuv->windowsize)/3;
  for (i = 0; !had_err && i < codons; i++) {
    double rvals[3] = {0.0,0.0,0.0};
    had_err = codon_usage_visitor_calc_window(cuv, gt_str_get(seq)+(i*3),
                                              rvals, err);
    acc[0] += rvals[0];
    acc[1] += rvals[1];
    acc[2] += rvals[2];
  }
  result[0] = acc[0]/codons;
  result[1] = acc[1]/codons;
  result[2] = acc[2]/codons;

  return had_err;
}

static int codon_usage_visitor_read_freqfile(GtCodonUsageVisitor *cuv,
                                             const char *filename,
                                             GtError *err)
{
  GtFile *file;
  GtStr *line;
  int had_err = 0, lineno = 1, codonno = 0;
  gt_assert(cuv && filename);

  file = gt_file_open(gt_file_mode_determine(filename), filename, "r", err);
  if (!file)
    had_err = -1;
  if (!had_err) {
    line = gt_str_new();
    while (!had_err && gt_str_read_next_line_generic(line, file) != EOF) {
      GT_UNUSED GtUword cnt, total;
      double logfreq;
      if (3 != sscanf(gt_str_get(line), GT_WU "\t" GT_WU "\t%lf",
                      &cnt, &total, &logfreq)) {
        had_err = -1;
        gt_error_set(err, "malformed entry in frequency file %s, line %d",
                     filename, lineno);
      }
      cuv->codon_freqs[codonno++] = logfreq;
      gt_str_reset(line);
      lineno++;
    }
    gt_str_delete(line);
    gt_file_delete(file);
  }
  return had_err;
}

static inline int codon_usage_visitor_determine_coding_frame(double *val)
{
  if (val[0] > val[1]) {
    if (val[0] > val[2]) {
      return 0;
    } else {
      return 2;
    }
  } else {
    if (val[1] > val[2]) {
      return 1;
    } else {
      return 2;
    }
  }
}

static int codon_usage_visitor_feature_node(GtNodeVisitor *nv,
                                            GtFeatureNode *fn,
                                            GtError *err)
{
  GtFeatureNodeIterator *fni;
  GtFeatureNode *node;
  GtCodonUsageVisitor *cuv;
  double val[3] = {0.0, 0.0, 0.0};
  cuv = codon_usage_visitor_cast(nv);
  int had_err = 0;

  fni = gt_feature_node_iterator_new(fn);
  while (!had_err && (node = gt_feature_node_iterator_next(fni))) {
    GtRange rng;
    char buf[BUFSIZ];
    int coding_frame, feature_frame;
    gt_str_reset(cuv->seq);

    /* CDS only */
    if (gt_feature_node_get_type(node) != gt_symbol("CDS"))
      continue;

    /* do not look at CDS shorter than window size */
    rng = gt_genome_node_get_range((GtGenomeNode*) node);
    if (gt_feature_node_get_strand(node) == GT_STRAND_REVERSE) {
      cuv->startpos = rng.start - 1;
    } else {
      cuv->startpos = rng.start - 1;
    }
    if (gt_range_length(&rng) <= cuv->windowsize)
      continue;

    if (gt_extract_feature_sequence(cuv->seq, (GtGenomeNode*) node,
                                    gt_symbol("CDS"),
                                    false, NULL, NULL, cuv->rm, err)) {
      had_err = -1;
    }
    if (!had_err)
      had_err = codon_usage_visitor_calc_value(cuv, cuv->seq, val, err);

    if (!had_err) {
      coding_frame = codon_usage_visitor_determine_coding_frame(val);
      feature_frame = ((int) (cuv->startpos) % 3);

      gt_log_log("pred: %d%c orig: %d%c: (%s) p0 %lf p1 %lf p2 %lf",
                 coding_frame,
                 GT_STRAND_CHARS[gt_feature_node_get_strand(node)],
                 feature_frame,
                 GT_STRAND_CHARS[gt_feature_node_get_strand(node)],
                 (coding_frame == feature_frame ? "COD" : "   "),
                 val[0], val[1], val[2]);

      snprintf(buf, BUFSIZ, "%d", ((int) (cuv->startpos) % 3));
      gt_feature_node_set_attribute(node, "frame", buf);
      snprintf(buf, BUFSIZ, "p0 %lf p1 %lf p2 %lf", val[0], val[1], val[2]);
      if (coding_frame == feature_frame
            && MAX3(val[0], val[1], val[2]) > cuv->coding_threshold) {
        gt_feature_node_set_attribute(node, "results", buf);
        gt_feature_node_set_score(node, MAX3(val[0], val[1], val[2]));
      }
    }
  }
  gt_feature_node_iterator_delete(fni);

  return had_err;
}

void codon_usage_visitor_free(GtNodeVisitor *nv)
{
  GtCodonUsageVisitor *cuv;
  if (!nv) return;
  cuv = codon_usage_visitor_cast(nv);
  gt_region_mapping_delete(cuv->rm);
  gt_str_delete(cuv->seq);
}

const GtNodeVisitorClass* gt_codon_usage_visitor_class()
{
  static GtNodeVisitorClass *nvc = NULL;
  gt_class_alloc_lock_enter();
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtCodonUsageVisitor),
                                    codon_usage_visitor_free,
                                    NULL,
                                    codon_usage_visitor_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}

GtNodeVisitor* gt_codon_usage_visitor_new(GtRegionMapping *rm,
                                          const char *freqtablefile,
                                          unsigned int windowsize,
                                          double coding_threshold,
                                          GtError *err)
{
  GtCodonUsageVisitor *cuv;
  GtNodeVisitor *nv;
  nv = gt_node_visitor_create(gt_codon_usage_visitor_class());
  cuv = codon_usage_visitor_cast(nv);

  cuv->rm = gt_region_mapping_ref(rm);
  cuv->seq = gt_str_new();
  cuv->windowsize = windowsize;
  cuv->coding_threshold = coding_threshold;

  if (codon_usage_visitor_read_freqfile(cuv, freqtablefile, err)) {
    gt_node_visitor_delete(nv);
    return NULL;
  }

  return nv;
}
