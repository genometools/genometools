/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GENOME_STREAM_REP_H
#define GENOME_STREAM_REP_H

#include <stdio.h>
#include "genome_stream.h"

struct GenomeStreamClass
{
  size_t size;
  GenomeNode* (*next_tree)(GenomeStream*, Log*);
  void         (*free)(GenomeStream*);
};

struct GenomeStream
{
  const GenomeStreamClass *c_class;
  GenomeNode *last_node;
  bool ensure_sorting;
};

GenomeStream* genome_stream_create(const GenomeStreamClass*,
                                    bool ensure_sorting);
void*          genome_stream_cast(const GenomeStreamClass*, GenomeStream*);

#endif
