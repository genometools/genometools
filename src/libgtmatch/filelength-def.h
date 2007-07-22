/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FILELENGTH_DEF_H
#define FILELENGTH_DEF_H
#include <inttypes.h>

typedef struct
{
  uint64_t length,
           effectivelength;
} Filelengthvalues;

#endif
