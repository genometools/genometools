/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GENOME_STREAM_H
#define GENOME_STREAM_H

#include <stdbool.h>
#include "error.h"
#include "genome_node.h"

/* the ``genome stream'' interface */
typedef struct GenomeStreamClass GenomeStreamClass;
typedef struct GenomeStream GenomeStream;

int   genome_stream_next_tree(GenomeStream*, GenomeNode**, Log*, Error*);
bool  genome_stream_is_sorted(GenomeStream*);
void  genome_stream_delete(GenomeStream*);

#endif
