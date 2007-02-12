/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef EXTRACTFEAT_VISITOR_H
#define EXTRACTFEAT_VISITOR_H

/* implements the ``genome visitor'' interface */
typedef struct Extractfeat_visitor Extractfeat_visitor;

#include <stdbool.h>
#include "genome_visitor.h"
#include "regionmapping.h"

const Genome_visitor_class* extractfeat_visitor_class(void);
Genome_visitor*             extractfeat_visitor_new_seqfile(Str *sequence_file,
                                                       Genome_feature_type type,
                                                       bool join,
                                                       bool translate);
/* takes ownership of the RegionMapping */
Genome_visitor*            extractfeat_visitor_new_regionmapping(RegionMapping*,
                                                       Genome_feature_type type,
                                                       bool join,
                                                       bool translate);

#endif
