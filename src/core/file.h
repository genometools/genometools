/*
  Copyright (c) 2026 Gordon Gremme <gordon@gremme.org>

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

#include <stdbool.h>
#include <stddef.h>
#include "core/file_api.h"

/* Internal block-wise read access to a <GtFile>, used to read lines without
   a function call (and a stream lock) per character.

   Only files which have been opened by <GtFile> itself via a path support
   this. <stdin> is excluded to keep its line-by-line behaviour for pipes and
   terminals, files created with <gt_file_new_from_fileptr()> are excluded
   because their <FILE*> can still be used directly by the caller.

   Bytes which are in the buffer are also returned by <gt_file_xfgetc()> and
   <gt_file_xread()>, so these functions can be mixed freely. */

/* Returns <true> if <file> supports the read buffer functions below. */
bool   gt_file_has_read_buffer(const GtFile *file);

/* Makes <*data> point to the unread bytes in the read buffer of <file>
   (refilling it, if it is empty) and returns their number. Returns 0 if the
   end of <file> has been reached. The bytes stay unread until
   <gt_file_read_buffer_consume()> is called, <*data> is valid until the next
   operation on <file>. */
size_t gt_file_read_buffer_peek(GtFile *file, const char **data);

/* Marks the first <nbytes> bytes returned by the last call to
   <gt_file_read_buffer_peek()> as read. */
void   gt_file_read_buffer_consume(GtFile *file, size_t nbytes);

#endif
