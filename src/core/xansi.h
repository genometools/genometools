/*
  Copyright (c) 2005-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "core/assert_api.h"

/*
  This module contains wrappers for the functions from the standard library we
  use. These functions always terminate the program if an error occurs.
  That is, one can use this functions without the need to check for errors.
*/

void   gt_xatexit(void (*function)(void));
void   gt_xfclose(FILE*);
void   gt_xfflush(FILE*);
int    gt_xfgetc(FILE*);
void   gt_xfgetpos(FILE*, fpos_t*);
FILE*  gt_xfopen(const char *path, const char *mode);
void   gt_xfputc(int, FILE*);
void   gt_xfputs(const char*, FILE*);
size_t gt_xfread(void *ptr, size_t size, size_t nmemb, FILE*);
void   gt_xfseek(FILE*, long offset, int whence);
void   gt_xfsetpos(FILE*, const fpos_t*);
void   gt_xfwrite(const void *ptr, size_t size, size_t nmemb, FILE*);
void   gt_xputchar(int);
void   gt_xputs(const char*);
void   gt_xremove(const char*);
char*  gt_xstrdup(const char*);
void   gt_xungetc(int, FILE*);
void   gt_xvfprintf(FILE *stream, const char *format, va_list ap);

#endif
