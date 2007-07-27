/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef SEQDESC_H
#define SEQDESC_H

#include "gqueue-def.h"

typedef struct
{
  Genericqueue *descptr;
  Arraychar headerbuffer;
} Sequencedescription;

#endif
