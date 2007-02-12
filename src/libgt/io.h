/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef IO_H
#define IO_H

#include <stdbool.h>
#include <stdio.h>

/* the I/O class */
typedef struct IO IO;

IO*           io_new(const char *path, const char *mode);
int           io_get_char(IO*, char*); /* returns -1 if no char is left */
void          io_unget_char(IO*, char);
bool          io_line_start(const IO*);
unsigned long io_get_line_number(const IO*);
const char*   io_get_filename(const IO*);
void          io_free(IO*);

#endif
