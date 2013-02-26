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
#include "core/fasta.h"
#include "core/ma.h"
#include "core/parseutils.h"
#include "core/range.h"
#include "core/str_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "extended/node_visitor_api.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/feature_type.h"
#include "ltr/ltrharvest_fasta_out_visitor.h"

struct GtLTRharvestFastaOutVisitor {
  const GtNodeVisitor parent_instance;
  GtFile* outfp;
  unsigned long width;
  bool inner;
  const GtEncseq *encseq;
};

#define gt_ltrharvest_fasta_out_visitor_cast(GV)\
        gt_node_visitor_cast(gt_ltrharvest_fasta_out_visitor_class(), GV)

static int gt_ltrharvest_fasta_out_visitor_feature_node(GtNodeVisitor *nv,
                                                        GtFeatureNode *fn,
                                                        GT_UNUSED GtError *err)
{
  GtLTRharvestFastaOutVisitor *lv;
  GtFeatureNodeIterator *fni;
  GtFeatureNode *curnode = NULL,
                *leftltr = NULL,
                *rightltr = NULL,
                *ltr_retrotrans = NULL;
  int had_err = 0,
      added_ltr = 0;
  GtRange rng,
          outrng;
  unsigned long seqnum = GT_UNDEF_ULONG;
  const char *fnt;
  lv = gt_ltrharvest_fasta_out_visitor_cast(nv);
  gt_assert(lv);
  gt_error_check(err);

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
  }
  gt_feature_node_iterator_delete(fni);

  /* check for invalid annotations */
  if (!had_err && (!leftltr || !rightltr)) {
    gt_error_set(err, "missing LTR annotations");
    had_err = -1;
  }

  /* determine sequence interval to output */
  if (lv->inner) {
    gt_assert(leftltr && rightltr);
    rng = gt_genome_node_get_range((GtGenomeNode*) leftltr);
    outrng.start = rng.end + 1;
    rng = gt_genome_node_get_range((GtGenomeNode*) rightltr);
    outrng.end = rng.start - 1;
  } else {
    gt_assert(ltr_retrotrans);
    rng = gt_genome_node_get_range((GtGenomeNode*) ltr_retrotrans);
    outrng = rng;
  }

  /* output FASTA sequences */
  if (!had_err && ltr_retrotrans != NULL) {
    if (outrng.start < outrng.end) {
      char *buf;
      const char *seqdesc;
      GtStr *desc;
      unsigned long startpos,
                    seqdesclen;
      gt_assert(seqnum != GT_UNDEF_ULONG
                  && seqnum < gt_encseq_num_of_sequences(lv->encseq));
      seqdesc = gt_encseq_description(lv->encseq, &seqdesclen, seqnum);
      desc = gt_str_new();
      gt_str_append_cstr_nt(desc, seqdesc, seqdesclen);
      gt_str_append_cstr(desc, " (dbseq-nr ");
      gt_str_append_ulong(desc, seqnum);
      gt_str_append_cstr(desc, ") [");
      gt_str_append_ulong(desc, outrng.start);
      gt_str_append_cstr(desc, ",");
      gt_str_append_ulong(desc, outrng.end);
      gt_str_append_cstr(desc, "]");
      buf = gt_calloc((size_t) gt_range_length(&outrng) + 1, sizeof (char));
      startpos = gt_encseq_seqstartpos(lv->encseq, seqnum);
      gt_encseq_extract_decoded(lv->encseq, buf,
                                startpos + outrng.start - 1,
                                startpos + outrng.end - 1);
      gt_fasta_show_entry(gt_str_get(desc), buf, gt_range_length(&outrng),
                          lv->width, lv->outfp);
      gt_free(buf);
      gt_str_delete(desc);
    } else {
      GtRange rootrng;
      rootrng = gt_genome_node_get_range((GtGenomeNode*) ltr_retrotrans);
      gt_warning("trying to output empty%s sequence for candidate at "
                 "%lu-%lu on sequence %lu",
                 (lv->inner ? " inner" : ""),
                 rootrng.start,
                 rootrng.end,
                 seqnum);
    }
  }

  return had_err;
}

const GtNodeVisitorClass* gt_ltrharvest_fasta_out_visitor_class(void)
{
  static const GtNodeVisitorClass *nvc = NULL;
  gt_class_alloc_lock_enter();
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtLTRharvestFastaOutVisitor),
                                   NULL,
                                   NULL,
                                   gt_ltrharvest_fasta_out_visitor_feature_node,
                                   NULL,
                                   NULL,
                                   NULL);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}

GtNodeVisitor* gt_ltrharvest_fasta_out_visitor_new(const GtEncseq *encseq,
                                                   bool inner,
                                                   unsigned long width,
                                                   GtFile *outfp)
{
  GtNodeVisitor *nv;
  GtLTRharvestFastaOutVisitor *lv;
  gt_assert(encseq && outfp);
  nv = gt_node_visitor_create(gt_ltrharvest_fasta_out_visitor_class());
  lv = gt_ltrharvest_fasta_out_visitor_cast(nv);
  gt_assert(lv);
  lv->inner = inner;
  lv->outfp = outfp;
  lv->width = width;
  lv->encseq = encseq;
  return nv;
}
