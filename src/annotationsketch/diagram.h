/*
  Copyright (c) 2007      Malte Mader <mmader@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007      Christin Schaerfer <cschaerfer@zbh.uni-hamburg.de>
  Copyright (c) 2007      Center for Bioinformatics, University of Hamburg

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

#ifndef DIAGRAM_H
#define DIAGRAM_H

#include "annotationsketch/diagram_api.h"
#include "core/hashmap.h"

typedef struct GtTracklineInfo {
  unsigned long total_lines,
                total_captionlines;
} GtTracklineInfo;

GtHashmap*    gt_diagram_get_tracks(const GtDiagram*);
void        gt_diagram_get_lineinfo(const GtDiagram*, GtTracklineInfo*);
int         gt_diagram_get_number_of_tracks(const GtDiagram*);
int         gt_diagram_unit_test(GtError*);

#endif
