/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/
/**
 * \file feature_visitor.h
 * \author Gordon Gremme <gremme@zbh.uni-hamburg.de>
 */

#ifndef FEATURE_VISITOR_H
#define FEATURE_VISITOR_H

/* implements the ``genome visitor'' interface */
typedef struct FeatureVisitor FeatureVisitor;

#include "libgtext/genome_visitor.h"
#include <libgtview/feature_index.h>

const GenomeVisitorClass* feature_visitor_class(void);
GenomeVisitor*            feature_visitor_new(FeatureIndex* , Env*);

#endif
