/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef XZLIB_H
#define XZLIB_H

#include <zlib.h>

/*
  This module contains wrappers for the functions from the zlib we use.
  These functions always terminate the program if an error occurs.
  That is, one can use this functions without the need to check for errors.
*/

gzFile xgzopen(const char *path, const char *mode);
void   xgzclose(gzFile);

#endif
