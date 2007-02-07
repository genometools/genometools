/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef ORF_H
#define ORF_H

#include "array.h"

/* the determined ORFs include the start and the stop codon */
void determine_ORFs(Array *ranges, unsigned int framenum,
                    const char *frame, unsigned long framelen);

#endif
