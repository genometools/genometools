/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FILEUTILS_H
#define FILEUTILS_H

#include <stdbool.h>
#include <stdio.h>
#include <libgtcore/str.h>

bool           file_exists(const char*);
/* returns 1 if the file with path a has a later modification time then the file
   with path b, 0 otherwise. */
bool           file_is_newer(const char *a, const char *b);
unsigned long  file_number_of_lines(FILE*);

/* set <path> to the dirname of <file>, if it has one, to "" otherwise */
void           file_dirname(Str *path, const char *file, Env*);

/* find 'file' in $PATH, if it has no dirname; set 'path' to dirname otherwise.
   sets 'path' to the empty string if 'file' could not be found in $PATH. */
int            file_find_in_path(Str *path, const char *file, Env*);

#endif
