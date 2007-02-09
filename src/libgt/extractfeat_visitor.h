/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef EXTRACTFEAT_VISITOR_H
#define EXTRACTFEAT_VISITOR_H

/* implements the ``genome visitor'' interface */
typedef struct Extractfeat_visitor Extractfeat_visitor;

#include "genome_visitor.h"

const Genome_visitor_class* extractfeat_visitor_class(void);
Genome_visitor*             extractfeat_visitor_new(Str *sequence_file,
                                                    Genome_feature_type type,
                                                    unsigned int join,
                                                    unsigned int translate);

#endif
