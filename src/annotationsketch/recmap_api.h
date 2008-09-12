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

#ifndef RECMAP_API_H
#define RECMAP_API_H

#include "extended/genome_feature.h"

typedef struct GtRecMap GtRecMap;

double                  gt_recmap_get_northwest_x(const GtRecMap*);
double                  gt_recmap_get_northwest_y(const GtRecMap*);
double                  gt_recmap_get_southeast_x(const GtRecMap*);
double                  gt_recmap_get_southeast_y(const GtRecMap*);
const GT_GenomeFeature* gt_recmap_get_genome_feature(const GtRecMap*);
bool                    gt_recmap_has_omitted_children(const GtRecMap*);

#endif
