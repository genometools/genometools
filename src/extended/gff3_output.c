/*
  Copyright (c) 2006, 2008 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006, 2008 Center for Bioinformatics, University of Hamburg

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

#include "core/assert_api.h"
#include "core/undef_api.h"
#include "extended/genome_node.h"
#include "extended/gff3_output.h"

void gt_gff3_output_leading(GtFeatureNode *fn, GtFile *outfp)
{
  GtGenomeNode *gn;
  gt_assert(fn);
  gn = (GtGenomeNode*) fn;
  gt_file_xprintf(outfp, "%s\t%s\t%s\t"GT_WU"\t"GT_WU"\t",
                     gt_str_get(gt_genome_node_get_seqid(gn)),
                     gt_feature_node_get_source(fn),
                     gt_feature_node_get_type(fn),
                     gt_genome_node_get_start(gn),
                     gt_genome_node_get_end(gn));
  if (gt_feature_node_score_is_defined(fn))
    gt_file_xprintf(outfp, "%.3g", gt_feature_node_get_score(fn));
  else
    gt_file_xfputc('.', outfp);
  gt_file_xprintf(outfp, "\t%c\t%c\t",
                     GT_STRAND_CHARS[gt_feature_node_get_strand(fn)],
                     GT_PHASE_CHARS[gt_feature_node_get_phase(fn)]);
}

void gt_gff3_output_leading_str(GtFeatureNode *fn, GtStr *outstr)
{
  GtGenomeNode *gn;
  gt_assert(fn && outstr);
  gn = (GtGenomeNode*) fn;
  gt_str_append_str(outstr, gt_genome_node_get_seqid(gn));
  gt_str_append_char(outstr, '\t');
  gt_str_append_cstr(outstr, gt_feature_node_get_source(fn));
  gt_str_append_char(outstr, '\t');
  gt_str_append_cstr(outstr, gt_feature_node_get_type(fn));
  gt_str_append_char(outstr, '\t');
  gt_str_append_ulong(outstr, gt_genome_node_get_start(gn));
  gt_str_append_char(outstr, '\t');
  gt_str_append_ulong(outstr, gt_genome_node_get_end(gn));
  gt_str_append_char(outstr, '\t');
  if (gt_feature_node_score_is_defined(fn)) {
    char buf[BUFSIZ];
    (void) snprintf(buf, BUFSIZ, "%.3g", gt_feature_node_get_score(fn));
    gt_str_append_cstr(outstr, buf);
  } else
    gt_str_append_char(outstr, '.');
  gt_str_append_char(outstr, '\t');
  gt_str_append_char(outstr, GT_STRAND_CHARS[gt_feature_node_get_strand(fn)]);
  gt_str_append_char(outstr, '\t');
  gt_str_append_char(outstr, GT_PHASE_CHARS[gt_feature_node_get_phase(fn)]);
  gt_str_append_char(outstr, '\t');
}