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

#ifndef FASTQ_H
#define FASTQ_H

#include "core/file.h"

#define GT_FASTQ_SEPARATOR_SEQ '@'
#define GT_FASTQ_SEPARATOR_QUAL '+'

/* Print a fastq entry with optional <description> and mandatory
   <sequence> and <qualities> to <outfp>. If <width> is != 0 the sequence and
   the qualities are formatted accordingly. If <repeat_description> is
   true, the description is output also before the qualities, otherwise
   only before the sequence. */
void gt_fastq_show_entry(const char *description, const char *sequence,
                         const char *qualities, unsigned long sequence_length,
                         unsigned long width, bool repeat_description,
                         GtFile *outfp);

#endif
