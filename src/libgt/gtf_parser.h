/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GTF_PARSER_H
#define GTF_PARSER_H

#include "env.h"
#include "queue.h"

/* This is a parser for gth GTF2.2 Gene Annotation Format as described at
   http://genes.cs.wustl.edu/GTF22.html

   It does not implement parsing of the following features given in the spec:

   - 5UTR
   - 3UTR
   - inter
   - inter_CNS
   - intron_CNS

*/

typedef struct GTF_parser GTF_parser;

GTF_parser* gtf_parser_new(void);
int         gtf_parser_parse(GTF_parser*, Queue *genome_nodes,
                             const char *filename, FILE*,
                             unsigned int be_tolerant, Env *err);
void        gtf_parser_delete(GTF_parser*);

#endif
