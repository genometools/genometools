/*
  Copyright (c) 2005-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2008 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
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

void   xatexit(void (*function)(void));
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
char*  xstrdup(const char*);
void   xungetc(int, FILE*);

#endif
