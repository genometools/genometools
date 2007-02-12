/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef EXTRACTFEAT_STREAM_H
#define EXTRACTFEAT_STREAM_H

#include <stdio.h>
#include "genome_stream.h"
#include "str.h"

/* implements the ``genome_stream'' interface */
typedef struct Extractfeat_stream Extractfeat_stream;

const Genome_stream_class* extractfeat_stream_class(void);

Genome_stream*             extractfeat_stream_new(Genome_stream*,
                                                  Str *sequence_file,
                                                  Genome_feature_type type,
                                                  bool join,
                                                  bool translate);

#endif
