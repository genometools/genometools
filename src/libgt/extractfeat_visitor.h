/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef EXTRACTFEAT_VISITOR_H
#define EXTRACTFEAT_VISITOR_H

/* implements the ``genome visitor'' interface */
typedef struct ExtractFeatVisitor ExtractFeatVisitor;

#include <stdbool.h>
#include "genome_visitor.h"
#include "regionmapping.h"

const GenomeVisitorClass* extractfeat_visitor_class(void);
GenomeVisitor*            extractfeat_visitor_new_seqfile(Str *sequence_file,
                                                          GenomeFeatureType,
                                                          bool join,
                                                          bool translate,
                                                          Env*);
/* takes ownership of the RegionMapping */
GenomeVisitor*            extractfeat_visitor_new_regionmapping(RegionMapping*,
                                                              GenomeFeatureType,
                                                                bool join,
                                                                bool translate,
                                                                Env*);

#endif
