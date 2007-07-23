/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef MINMAX_H
#define MINMAX_H

#ifndef MAX
#define MAX(X,Y) (((X) > (Y)) ? (X) : (Y))
#endif

#ifndef MIN
#define MIN(X,Y) (((X) < (Y)) ? (X) : (Y))
#endif

#ifndef MIN3
#define MIN3(a, b, c) (((a)<(b))?((a)<(c)?(a):(c)):((b)<(c)?(b):(c)))
#endif

#endif
