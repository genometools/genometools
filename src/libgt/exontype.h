/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef EXONTYPE_H
#define EXONTYPE_H

typedef enum {
  EXONTYPE_SINGLE,
  EXONTYPE_INITIAL,
  EXONTYPE_INTERNAL,
  EXONTYPE_TERMINAL,
  EXONTYPE_UNDETERMINED
} ExonType;

#endif
