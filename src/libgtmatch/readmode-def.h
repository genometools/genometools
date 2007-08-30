/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef READMODE_DEF_H
#define READMODE_DEF_H

typedef enum
{
  Forwardmode = 0,
  Reversemode,
  Complementmode,
  Reversecomplementmode
} Readmode;

#define ISDIRREVERSE(R) ((R) == Reversemode || (R) == Reversecomplementmode)

#endif
