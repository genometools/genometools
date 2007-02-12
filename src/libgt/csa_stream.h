/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef CSA_STREAM_H
#define CSA_STREAM_H

#include <stdio.h>
#include "genome_stream.h"

/* implements the ``genome_stream'' interface */
typedef struct Csa_stream Csa_stream;

const GenomeStreamClass* csa_stream_class(void);
GenomeStream*             csa_stream_new(GenomeStream*,
                                          unsigned long join_length);

#endif
