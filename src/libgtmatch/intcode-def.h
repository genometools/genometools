/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef INTCODE_DEF_H
#define INTCODE_DEF_H
#include <stdbool.h>

#define PREFIXLENBITS   4

typedef unsigned int Codetype;      /* \Typedef{Codetype} */

typedef struct
{
  unsigned int specialpos;
  bool defined;
} Firstspecialpos;

#endif
