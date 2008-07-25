/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include "libgtcore/ma.h"
#include "libgtext/genome_feature.h"
#include "libgtext/genome_node.h"
#include "libgtview/recmap.h"

RecMap* recmap_create(double nw_x, double nw_y,
                      double se_x, double se_y,
                      GenomeFeature *gf)
{
  RecMap *rm;
  rm = ma_calloc(1, sizeof (RecMap));
  assert(rm);
  rm->nw_x = nw_x;
  rm->nw_y = nw_y;
  rm->se_x = se_x;
  rm->se_y = se_y;
  rm->gf = (GenomeFeature*) genome_node_ref((GenomeNode*) gf);
  return rm;
}

int recmap_format_html_imagemap_coords(RecMap *rm, char *buf, size_t n)
{
  assert(rm && buf);
  return snprintf(buf, n, "%.0f,%.0f,%.0f,%.0f", rm->nw_x, rm->nw_y,
                                                 rm->se_x, rm->se_y);
}

void recmap_delete(RecMap *rm)
{
  if (!rm) return;
  genome_node_delete((GenomeNode*) rm->gf);
  ma_free(rm);
}
