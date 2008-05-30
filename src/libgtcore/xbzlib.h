/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef XBZLIB_H
#define XBZLIB_H

#include <bzlib.h>

/*
  This module contains wrappers for the functions from the bz2lib we use.
  These functions always terminate the program if an error occurs.
  That is, one can use this functions without the need to check for errors.
*/

BZFILE* xbzopen(const char *path, const char *mode);
/* Returns next character from <bzfile> or EOF, if end-of-file is reached. */
int     xbzfgetc(BZFILE *bzfile);
void    xbzfputc(int, BZFILE*);
void    xbzfputs(const char*, BZFILE*);
/* Returns num of read bytes. */
int     xbzread(BZFILE*, void *buf, unsigned len);
void    xbzwrite(BZFILE*, void *buf, unsigned len);
void    xbzrewind(BZFILE**, const char *orig_path, const char *orig_mode);

#endif
