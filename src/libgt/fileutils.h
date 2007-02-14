/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FILEUTILS_H
#define FILEUTILS_H

#include <stdbool.h>
#include <stdio.h>
#include "str.h"

bool          file_exists(const char*);
/* returns 1 if the file with path a has a later modification time then the file
   with path b, 0 otherwise. */
bool          file_is_newer(const char *a, const char *b);
unsigned long file_number_of_lines(FILE*);

/* return dirname of 'file', if it has one, NULL otherwise */
Str*          file_dirname(const char *file);

/* find 'file' in $PATH, if it has no dirname; return dirname otherwise.
   returns NULL is 'file' could not be found in $PATH. */
Str*          file_find_in_path(const char *file);

#endif
