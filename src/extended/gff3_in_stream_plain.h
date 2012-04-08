/*
  Copyright (c) 2006-2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef GFF3_IN_STREAM_PLAIN_H
#define GFF3_IN_STREAM_PLAIN_H

#include <stdio.h>
#include "extended/gff3_in_stream_plain.h"
#include "extended/node_stream_api.h"
#include "extended/type_checker_api.h"

/* Implements the <GtNodeStream> interface. */
typedef struct GtGFF3InStreamPlain GtGFF3InStreamPlain;

const GtNodeStreamClass* gt_gff3_in_stream_plain_class(void);

GtNodeStream* gt_gff3_in_stream_plain_new_unsorted(int num_of_files,
                                                   const char **filenames);
GtNodeStream* gt_gff3_in_stream_plain_new_sorted(const char *filename);
void          gt_gff3_in_stream_plain_check_id_attributes(GtGFF3InStreamPlain*);
void          gt_gff3_in_stream_plain_check_region_boundaries(
                                                          GtGFF3InStreamPlain*);
void          gt_gff3_in_stream_plain_do_not_check_region_boundaries(
                                                          GtGFF3InStreamPlain*);
void          gt_gff3_in_stream_plain_enable_tidy_mode(GtNodeStream*);
void          gt_gff3_in_stream_plain_enable_strict_mode(GtNodeStream*);
void          gt_gff3_in_stream_plain_show_progress_bar(GtGFF3InStreamPlain*);
void          gt_gff3_in_stream_plain_set_type_checker(GtNodeStream*,
                                                       GtTypeChecker*);
GtStrArray*   gt_gff3_in_stream_plain_get_used_types(GtNodeStream*);
void          gt_gff3_in_stream_plain_set_offset(GtNodeStream*, long);
int           gt_gff3_in_stream_plain_set_offsetfile(GtNodeStream*, GtStr*,
                                                     GtError*);

#endif
