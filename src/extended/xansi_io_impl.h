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

#ifndef XANSI_IO_IMPL_H
#define XANSI_IO_IMPL_H

#include <stdio.h>
#include <stdlib.h>

#include "core/assert_api.h"
#include "core/unused_api.h"

GT_UNUSED static inline void gt_xansi_io_xfwrite(void *ptr,
                                       size_t size,
                                       size_t nmemb,
                                       FILE *stream)
{
  if (nmemb != fwrite((const void*) ptr, size, nmemb, stream)) {
    perror("gt_xansi_io_xfwrite failed to write to file");
    exit(EXIT_FAILURE);
  }
}

GT_UNUSED static inline void gt_xansi_io_xfread(void *ptr,
                                      size_t size,
                                      size_t nmemb,
                                      FILE *stream)
{
  if (nmemb != fread(ptr, size, nmemb, stream)) {
    gt_assert(feof(stream) == 0);
    if (ferror(stream) != 0)
      perror("gt_xansi_io_xfread failed to read from file");
    exit(EXIT_FAILURE);
  };
}
#endif
