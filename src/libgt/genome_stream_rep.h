/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GENOME_STREAM_REP_H
#define GENOME_STREAM_REP_H

#include <stdio.h>
#include "genome_stream.h"

struct Genome_stream_class
{
  size_t size;
  GenomeNode* (*next_tree)(Genome_stream*, Log*);
  void         (*free)(Genome_stream*);
};

struct Genome_stream
{
  const Genome_stream_class *c_class;
  GenomeNode *last_node;
  bool ensure_sorting;
};

Genome_stream* genome_stream_create(const Genome_stream_class*,
                                    bool ensure_sorting);
void*          genome_stream_cast(const Genome_stream_class*, Genome_stream*);

#endif
