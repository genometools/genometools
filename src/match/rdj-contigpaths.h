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

#ifndef RDJ_CONTIGPATHS_H
#define RDJ_CONTIGPATHS_H

#include "core/file.h"
#include "core/encseq.h"
#include "core/logger.h"
#include "core/error_api.h"

typedef uint32_t GtContigpathElem;
#define GT_CONTIGPATH_ELEM_MAX (GtContigpathElem)UINT32_MAX

int gt_contigpaths_to_fasta(const char *indexname,
    const char *contigpaths_suffix, const char *fasta_suffix,
    const GtEncseq *encseq, unsigned long min_contig_length, bool showpaths,
    size_t buffersize, GtLogger *logger, GtError *err);

#endif
