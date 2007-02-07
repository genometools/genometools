/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FILEUTILS_H
#define FILEUTILS_H

#include <stdio.h>

unsigned int  file_exists(const char*);
/* returns 1 if the file with path a has a later modification time then the file
   with path b, 0 otherwise. */
unsigned int  file_is_newer(const char *a, const char *b);
unsigned long file_number_of_lines(FILE*);

#endif
