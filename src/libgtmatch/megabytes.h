/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef MEGABYTES_H
#define MEGABYTES_H

/*
  The following macro transforms bytes into megabytes.
*/

#define MEGABYTES(V)  ((double) (V)/((((unsigned long) 1) << 20) - 1))

#endif
