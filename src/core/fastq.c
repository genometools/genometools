/*
  Copyright (c) 2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "core/fastq.h"
#include "core/xansi_api.h"

static void gt_fastq_show_buffer(char separator, const char *description,
    const char *buffer, unsigned long buffer_length, unsigned long width,
    GtFile *outfp)
{
  unsigned long i, current_length;
  gt_file_xfputc(separator, outfp);
  if (description != NULL)
    gt_file_xfputs(description, outfp);
  gt_file_xfputc('\n', outfp);
  for (i = 0, current_length = 0; i < buffer_length;
       i++, current_length++) {
    if (width && current_length == width) {
      gt_file_xfputc('\n', outfp);
      current_length = 0;
    }
    gt_file_xfputc(buffer[i], outfp);
  }
  gt_file_xfputc('\n', outfp);
}

void gt_fastq_show_entry(const char *description, const char *sequence,
                         const char *qualities, unsigned long sequence_length,
                         unsigned long width, bool repeat_description,
                         GtFile *outfp)
{
  gt_assert(sequence);
  gt_assert(qualities);
  gt_fastq_show_buffer(GT_FASTQ_SEPARATOR_SEQ, description, sequence,
      sequence_length, width, outfp);
  gt_fastq_show_buffer(GT_FASTQ_SEPARATOR_QUAL,
      repeat_description ? description : NULL, qualities,
      sequence_length, width, outfp);
}
