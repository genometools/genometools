/*
  Copyright (c) 2006-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/cstr_table_api.h"
#include "core/file_api.h"
#include "core/queue_api.h"
#include "core/range_api.h"
#include "core/strand_api.h"
#include "extended/type_checker_api.h"

#define GT_GFF_VERSION         3
#define GT_GFF_VERSION_PREFIX  "##gff-version"
#define GT_GFF_FASTA_DIRECTIVE "##FASTA"
#define GT_GFF_SEQUENCE_REGION "##sequence-region"
#define GT_GFF_TERMINATOR      "###"

#define GT_GFF_ID              "ID"
#define GT_GFF_PARENT          "Parent"
#define GT_GFF_TARGET          "Target"

typedef struct GtGFF3Parser GtGFF3Parser;

GtGFF3Parser* gt_gff3_parser_new(GtTypeChecker*);
void          gt_gff3_parser_check_id_attributes(GtGFF3Parser*);
void          gt_gff3_parser_set_offset(GtGFF3Parser*, long);
int           gt_gff3_parser_set_offsetfile(GtGFF3Parser*, GtStr*, GtError*);
void          gt_gff3_parser_enable_tidy_mode(GtGFF3Parser*);
int           gt_gff3_parser_parse_target_attributes(const char *values,
                                                     unsigned long
                                                     *num_of_targets,
                                                     GtStr *first_target_id,
                                                     GtRange
                                                     *first_target_range,
                                                     GtStrand
                                                     *first_target_strand,
                                                     const char *filename,
                                                     unsigned int line_number,
                                                     GtError*);
void          gt_gff3_parser_parse_all_target_attributes(const char *values,
                                                         GtStrArray *target_ids,
                                                         GtArray *target_ranges,
                                                         GtArray
                                                         *target_strands);
int           gt_gff3_parser_parse_genome_nodes(GtGFF3Parser*,
                                                int *status_code,
                                                GtQueue *genome_nodes,
                                                GtCstrTable *used_types,
                                                GtStr *filenamestr,
                                                unsigned long long *line_number,
                                                GtFile *fpin,
                                                GtError *err);
/* Reset the GFF3 parser (necessary if the processed input file is switched). */
void          gt_gff3_parser_reset(GtGFF3Parser*);
void          gt_gff3_parser_delete(GtGFF3Parser*);

#endif
