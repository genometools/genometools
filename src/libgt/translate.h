/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef TRANSLATE_H
#define TRANSLATE_H

#include "str.h"

void translate_dna(Str*, const char*, unsigned long, unsigned int frame);

#endif
