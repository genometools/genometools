/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GENOME_NODE_REP_H
#define GENOME_NODE_REP_H

#include <stdio.h>
#include "dlist.h"
#include "genome_node.h"

/* the ``genome node'' interface */
struct Genome_node_class
{
  size_t size;
  void  (*free)(Genome_node*);
  Str*  (*get_seqid)(Genome_node*);
  Str*  (*get_idstr)(Genome_node*);
  Range (*get_range)(Genome_node*);
  void  (*set_range)(Genome_node*, Range);
  void  (*set_seqid)(Genome_node*, Str*);
  void  (*set_source)(Genome_node*, Str*);
  void  (*set_phase)(Genome_node*, Phase);
  void  (*accept)(Genome_node*, Genome_visitor*, Log*);
};

struct Genome_node
{
  const Genome_node_class *c_class;
  const char *filename;
  unsigned long line_number;
  Dlist *children;
  unsigned int reference_count;
};

void         genome_node_class_init(Genome_node_class*, size_t, ...);
Genome_node* genome_node_create(const Genome_node_class*, const char *filename,
                                unsigned long line_number);

#endif
