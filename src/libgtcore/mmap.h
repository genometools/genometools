/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef MMAP_H
#define MMAP_H

#include <stdio.h>

/* This module contains hight level mmap(2) interfaces. */

void* mmap_read(const char *path, size_t *len);
void* mmap_write(const char *path, size_t *len);
void* xmmap_read(const char *path, size_t *len);
void* xmmap_write(const char *path, size_t *len);

#endif
