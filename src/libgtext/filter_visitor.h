/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef FILTER_VISITOR_H
#define FILTER_VISITOR_H

/* implements the ``genome visitor'' interface */
typedef struct FilterVisitor FilterVisitor;

#include "libgtext/genome_visitor.h"

const GenomeVisitorClass* filter_visitor_class(void);
/* If <strand> is != NUM_OF_STRAND_TYPES, then each genome feature must have
   strand <strand>. */
GenomeVisitor*            filter_visitor_new(Str *seqid, Str *typefilter,
                                             Range contain_range,
                                             Range overlap_range, Strand strand,
                                             Strand targetstrand, bool has_CDS,
                                             unsigned long max_gene_length,
                                             unsigned long max_gene_num,
                                             double min_gene_score,
                                             double
                                             min_average_splice_site_prob);
unsigned long             filter_visitor_node_buffer_size(GenomeVisitor*);
GenomeNode*               filter_visitor_get_node(GenomeVisitor*);

#endif
