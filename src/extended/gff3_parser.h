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

#ifndef GFF3_PARSER_H
#define GFF3_PARSER_H

#include "extended/gff3_parser_api.h"

void gt_gff3_parser_enable_strict_mode(GtGFF3Parser*);
int  gt_gff3_parser_set_offsetfile(GtGFF3Parser*, GtStr*, GtError*);
int  gt_gff3_parser_parse_target_attributes(const char *values,
                                            unsigned long *num_of_targets,
                                            GtStr *first_target_id,
                                            GtRange *first_target_range,
                                            GtStrand *first_target_strand,
                                            const char *filename,
                                            unsigned int line_number, GtError*);
int  gt_gff3_parser_parse_all_target_attributes(const char *values, bool tidy,
                                                GtStrArray *target_ids,
                                                GtArray *target_ranges,
                                                GtArray *target_strands,
                                                const char *filename,
                                                unsigned int line_number,
                                                GtError *err);
void gt_gff3_parser_build_target_str(GtStr *target, GtStrArray *target_ids,
                                     GtArray *target_ranges,
                                     GtArray *target_strands);

#endif
