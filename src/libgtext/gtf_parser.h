/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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

#ifndef GTF_PARSER_H
#define GTF_PARSER_H

#include "libgtcore/queue.h"

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
                             Str *filenamestr, FILE*, unsigned int be_tolerant,
                             Error*);
void        gtf_parser_delete(GTF_parser*);

#endif
