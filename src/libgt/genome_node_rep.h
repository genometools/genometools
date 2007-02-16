/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GENOME_NODE_REP_H
#define GENOME_NODE_REP_H

#include <stdio.h>
#include "dlist.h"
#include "genome_node.h"

/* the ``genome node'' interface */
struct GenomeNodeClass
{
  size_t size;
  void  (*free)(GenomeNode*);
  Str*  (*get_seqid)(GenomeNode*);
  Str*  (*get_idstr)(GenomeNode*);
  Range (*get_range)(GenomeNode*);
  void  (*set_range)(GenomeNode*, Range);
  void  (*set_seqid)(GenomeNode*, Str*);
  void  (*set_source)(GenomeNode*, Str*);
  void  (*set_phase)(GenomeNode*, Phase);
  int   (*accept)(GenomeNode*, GenomeVisitor*, Log*, Error*);
};

struct GenomeNode
{
  const GenomeNodeClass *c_class;
  const char *filename;
  unsigned long line_number;
  Dlist *children;
  unsigned int reference_count;
};

void        genome_node_class_init(GenomeNodeClass*, size_t, ...);
GenomeNode* genome_node_create(const GenomeNodeClass*, const char *filename,
                               unsigned long line_number);

#endif
