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
#include "libgtcore/undef.h"
#include "libgtext/gff3_output.h"

void gff3_output_leading(GenomeFeature *gf, GenFile *outfp)
{
  GenomeNode *gn;
  GenomeFeatureType type;
  double score;

  assert(gf);

  gn = (GenomeNode*) gf;
  type = genome_feature_get_type(gf);

  genfile_xprintf(outfp, "%s\t%s\t%s\t%lu\t%lu\t",
                  str_get(genome_node_get_seqid(gn)),
                  genome_feature_get_source(gf),
                  genome_feature_type_get_cstr(type),
                  genome_node_get_start(gn),
                  genome_node_get_end(gn));
  score = genome_feature_get_score(gf);
  if (score == UNDEF_DOUBLE)
    genfile_xfputc('.', outfp);
  else
    genfile_xprintf(outfp, "%.3f", score);
  genfile_xprintf(outfp, "\t%c\t%c\t",
                  STRANDCHARS[genome_feature_get_strand(gf)],
                  PHASECHARS[genome_feature_get_phase(gf)]);
}
