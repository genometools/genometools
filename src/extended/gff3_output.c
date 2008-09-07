/*
  Copyright (c) 2006, 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <assert.h>
#include "core/undef.h"
#include "extended/gff3_output.h"

void gff3_output_leading(GenomeFeature *gf, GenFile *outfp)
{
  GT_GenomeNode *gn;
  GenomeFeatureType *type;

  assert(gf);

  gn = (GT_GenomeNode*) gf;
  type = genome_feature_get_type(gf);

  genfile_xprintf(outfp, "%s\t%s\t%s\t%lu\t%lu\t",
                  str_get(gt_genome_node_get_seqid(gn)),
                  genome_feature_get_source(gf),
                  genome_feature_type_get_cstr(type),
                  gt_genome_node_get_start(gn),
                  gt_genome_node_get_end(gn));
  if (genome_feature_score_is_defined(gf))
    genfile_xprintf(outfp, "%.3g", genome_feature_get_score(gf));
  else
    genfile_xfputc('.', outfp);
  genfile_xprintf(outfp, "\t%c\t%c\t",
                  STRANDCHARS[genome_feature_get_strand(gf)],
                  PHASECHARS[genome_feature_get_phase(gf)]);
}
