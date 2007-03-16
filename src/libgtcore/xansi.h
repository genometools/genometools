/*
  Copyright (c) 2005-2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef XANSI_H
#define XANSI_H

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
  This module contains wrappers for the functions from the standard library we
  use. These functions always terminate the program if an error occurs.
  That is, one can use this functions without the need to check for errors.
*/

void*  xcalloc(size_t nmemb, size_t size);
void   xfclose(FILE*);
void   xfflush(FILE*);
int    xfgetc(FILE*);
void   xfgetpos(FILE*, fpos_t*);
FILE*  xfopen(const char *path, const char *mode);
void   xfputc(int, FILE*);
void   xfputs(const char*, FILE*);
size_t xfread(void *ptr, size_t size, size_t nmemb, FILE*);
void   xfseek(FILE*, long offset, int whence);
void   xfsetpos(FILE*, const fpos_t*);
void   xfwrite(const void *ptr, size_t size, size_t nmemb, FILE*);
void*  xmalloc(size_t size);
void   xputchar(int);
void   xputs(const char*);
void*  xrealloc(void *ptr, size_t size);
void   xremove(const char*);
char*  xtmpnam(char*);
void   xungetc(int, FILE*);

#endif
