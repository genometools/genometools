/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GFF3_PARSER_H
#define GFF3_PARSER_H

#include "queue.h"

#define GFF_VERSION_STRING  "##gff-version   3"
#define GFF_SEQUENCE_REGION "##sequence-region"
#define GFF_TERMINATOR      "###"

#define ID_STRING           "ID"
#define PARENT_STRING       "Parent"

typedef struct GFF3Parser GFF3Parser;

GFF3Parser* gff3_new(void);
void        gff3_set_offset(GFF3Parser*, long);
int         gff3_parse_genome_nodes(GFF3Parser*,
                                    Queue *genome_nodes,
                                    const char *filename,
                                    unsigned long *line_number,
                                    FILE *fpin);
/* resets the GFF3 parser (necessary if the processed input file is switched) */
void         gff3_reset(GFF3Parser*);
void         gff3_free(GFF3Parser*);

#endif
