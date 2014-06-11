/*
  Copyright (c) 2014 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2014 Center for Bioinformatics, University of Hamburg

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

#ifndef IO_FUNCTION_POINTERS_H
#define IO_FUNCTION_POINTERS_H

#include <stdio.h>
#include <stdlib.h>

#include "core/error_api.h"

/* IO function, either reads or writes to file <stream>. Either one of the two
   functions below. */
typedef int (*GtIOFunc)(void *ptr, size_t size, size_t nmemb, FILE *stream,
                        GtError *err);

/* Wrapper around fwrite(). Returns value other than 0 on error, and sets <err>
   accordingly. */
int gt_io_error_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream,
                       GtError *err);

#define gt_io_error_fwrite_one(element, fp, err) \
  gt_io_error_fwrite(&element, sizeof (element), (size_t) 1, fp, err)

/* Wrapper around fread(). Returns value other than 0 on error, and sets <err>
   accordingly. Less elements read than specified by <nmemb> is assumed to be
   always an error, even when reaching end of file. */
int gt_io_error_fread(void *ptr, size_t size, size_t nmemb, FILE *stream,
                      GtError *err);

#define gt_io_error_fread_one(element, fp, err) \
  gt_io_error_fread(&element, sizeof (element), (size_t) 1, fp, err)

#endif
