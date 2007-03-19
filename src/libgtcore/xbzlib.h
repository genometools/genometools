/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef XBZLIB_H
#define XBZLIB_H

#include <bzlib.h>

/*
  This module contains wrappers for the functions from the bz2lib we use.
  These functions always terminate the program if an error occurs.
  That is, one can use this functions without the need to check for errors.
*/

BZFILE* xbzopen(const char *path, const char *mode);
        /* returns num of read bytes */
int     xbzread(BZFILE*, void *buf, unsigned len);
void    xbzrewind(BZFILE**, const char *orig_path, const char *orig_mode);

#endif
