/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef EXTRACTFEAT_STREAM_H
#define EXTRACTFEAT_STREAM_H

#include <stdio.h>
#include "genome_stream.h"
#include "regionmapping.h"
#include "str.h"

/* implements the ``genome_stream'' interface */
typedef struct Extractfeat_stream Extractfeat_stream;

const Genome_stream_class* extractfeat_stream_class(void);

/* create a plain Extractfeat_stream */
Genome_stream*             extractfeat_stream_new(Genome_stream*,
                                                  Genome_feature_type type,
                                                  bool join,
                                                  bool translate);

/* exactly one of the following two functions has to be called to make an
   Extractfeat_stream usable */

/* use the file named ``seqfile'' as sequence source */
void                       extractfeat_stream_use_sequence_file(Genome_stream*,
                                                                Str *seqfile);

/* use the given RegionMapping (takes ownership) as sequence file source */
void                       extractfeat_stream_use_region_mapping(Genome_stream*,
                                                                RegionMapping*);

#endif
