/*
  Copyright (c) 2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#ifndef RDJ_CONTIGS_WRITER_H
#define RDJ_CONTIGS_WRITER_H

#include "core/encseq.h"
#include "core/file.h"
#include "core/logger.h"

typedef struct GtContigsWriter GtContigsWriter;

GtContigsWriter *gt_contigs_writer_new(const GtEncseq *reads, GtFile *outfp,
    bool showpaths);
void gt_contigs_writer_delete(GtContigsWriter *contigs_writer);

void gt_contigs_writer_start(GtContigsWriter *contigs_writer,
    unsigned long seqnum);

void gt_contigs_writer_append(GtContigsWriter *contigs_writer,
    unsigned long seqnum, unsigned long nofchars);

void gt_contigs_writer_write(GtContigsWriter *contigs_writer);

void gt_contigs_writer_abort(GtContigsWriter *contigs_writer);

void gt_contigs_writer_show_stats(GtContigsWriter *contigs_writer,
    GtLogger *logger);

#endif
