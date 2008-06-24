/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "libgtcore/queue.h"

#define GFF_VERSION         3
#define GFF_VERSION_PREFIX  "##gff-version"
#define GFF_SEQUENCE_REGION "##sequence-region"
#define GFF_TERMINATOR      "###"

#define ID_STRING           "ID"
#define PARENT_STRING       "Parent"

typedef struct GFF3Parser GFF3Parser;

GFF3Parser* gff3parser_new(bool checkids);
void        gff3parser_set_offset(GFF3Parser*, long);
int         gff3parser_set_offsetfile(GFF3Parser*, Str*, Error*);
int         gff3parser_parse_target_attributes(const char *values,
                                               unsigned long *num_of_targets,
                                               Str *first_target_id,
                                               Range *first_target_range,
                                               Strand *first_target_strand,
                                               const char *filename,
                                               unsigned long line_number,
                                               Error*);
int         gff3parser_parse_genome_nodes(int *status_code, GFF3Parser*,
                                          Queue *genome_nodes,
                                          Str *filenamestr,
                                          unsigned long long *line_number,
                                          GenFile *fpin, Error*);
/* resets the GFF3 parser (necessary if the processed input file is switched) */
void         gff3parser_reset(GFF3Parser*);
void         gff3parser_delete(GFF3Parser*);

#endif
