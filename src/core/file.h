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

#ifndef FILE_H
#define FILE_H

#include <stdlib.h>
#include "core/file_api.h"

typedef enum {
  GT_FILE_MODE_UNCOMPRESSED,
  GT_FILE_MODE_GZIP,
  GT_FILE_MODE_BZIP2
} GtFileMode;

/* Returns <GT_FILE_MODE_GZIP> if file with <path> ends with '.gz',
   <GT_FILE_MODE_BZIP2> if it ends with '.bz2', and <GT_FILE_MODE_UNCOMPRESSED>
   otherwise. */
GtFileMode  gt_file_mode_determine(const char *path);

/* Returns ".gz" if <mode> is GFM_GZIP, ".bz2" if <mode> is GFM_BZIP2, and ""
   otherwise. */
const char* gt_file_mode_suffix(GtFileMode mode);

/* Returns the length of the ``basename'' of <path>. That is, the length of path
   without '.gz' or '.bz2' suffixes. */
size_t      gt_file_basename_length(const char *path);

/* Create a new GtFile object and open the underlying file handle, returns
   NULL and sets <err> if the file <path> could not be opened. */
GtFile*     gt_file_open(GtFileMode, const char *path, const char *mode,
                         GtError*);

/* Create a new GtFile object and open the underlying file handle, abort if
   the file <path> does not exist. The <file_mode> has to be given
   explicitly. */
GtFile*     gt_file_xopen_file_mode(GtFileMode file_mode, const char *path,
                                    const char *mode);

/* Create a new GtFile object and open the underlying file handle. Aborts if
   the file <path> could not be opened. The GtFileMode is determined
   automatically via gt_file_mode_determine(path). */
GtFile*     gt_file_xopen(const char *path, const char *mode);

/* Returns the mode of the given <file>. */
GtFileMode  gt_file_mode(const GtFile *file);

/* Unget character <c> to <file> (which obviously cannot be <NULL>).
   Can only be used once at a time. */
void        gt_file_unget_char(GtFile *file, char c);

#endif
