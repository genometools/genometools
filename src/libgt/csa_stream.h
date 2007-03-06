/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef CSA_STREAM_H
#define CSA_STREAM_H

#include <stdio.h>
#include <libgt/genome_stream.h>

/* implements the ``genome_stream'' interface */
typedef struct CSAStream CSAStream;

const GenomeStreamClass* csa_stream_class(void);
GenomeStream*            csa_stream_new(GenomeStream*,
                                        unsigned long join_length, Env*);

#endif
