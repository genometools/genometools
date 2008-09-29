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

#ifndef XZLIB_H
#define XZLIB_H

#include <zlib.h>

/*
  This module contains wrappers for the functions from the zlib we use.
  These functions always terminate the program if an error occurs.
  That is, one can use this functions without the need to check for errors.
*/

gzFile gt_xgzopen(const char *path, const char *mode);
/* Returns next character from <file> or EOF, if end-of-file is reached. */
int    gt_xgzfgetc(gzFile file);
void   gt_xgzfputc(int, gzFile);
void   gt_xgzfputs(const char*, gzFile);
/* Returns num of read bytes. */
int    gt_xgzread(gzFile, void *buf, unsigned len);
void   gt_xgzwrite(gzFile, void *buf, unsigned len);
void   gt_xgzrewind(gzFile);
void   gt_xgzclose(gzFile);

#endif
