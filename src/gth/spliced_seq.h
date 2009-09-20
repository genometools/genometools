/*
  Copyright (c) 2004-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2004-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef SPLICED_SEQ_H
#define SPLICED_SEQ_H

#include "core/array_api.h"
#include "core/file.h"

typedef struct {
  GtArray *ranges;                /* the ranges of the genomic sequence which
                                     have been used for processing */
  const unsigned char *origseq;   /* the original sequence */
        unsigned char *splicedseq;/* the processed sequence */
  unsigned long splicedseqlen,    /* the length of the processed sequence */
                *positionmapping; /* maps the positions on the spliced sequence
                                     back to the original sequence */
} GthSplicedSeq;

GthSplicedSeq* gth_spliced_seq_new(const unsigned char *sequence,
                                   GtArray *ranges);
GthSplicedSeq* gth_spliced_seq_new_with_comments(const unsigned char *sequence,
                                                 GtArray *ranges, bool comments,
                                                 GtFile *outfp);
void           gth_spliced_seq_delete(GthSplicedSeq*);
bool           gth_spliced_seq_pos_is_border(const GthSplicedSeq*,
                                             unsigned long position);
unsigned long  gth_spliced_seq_border_length(const GthSplicedSeq*,
                                             unsigned long position);
unsigned long  gth_spliced_seq_num_of_borders(const GthSplicedSeq*);
/* The following function computes the spliced sequence position of a given
   original position <origpos> of a position mapping <positionmapping> of length
   <splicedseqlen> via binary search. */
unsigned long  gth_spliced_seq_orig_to_spliced_pos(const GthSplicedSeq*,
                                                   unsigned long orig_pos);

#endif
