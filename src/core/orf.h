/*
  Copyright (c) 2006, 2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006       Center for Bioinformatics, University of Hamburg

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

#ifndef ORF_H
#define ORF_H

#include "core/range_api.h"

typedef void (*GtORFProcessor)(void *data, GtRange *orf, unsigned long framenum,
                               const char *frame, bool ends_with_stop_codon);

/* Determine all ORFs in the given <frame> of length <framelength> and frame
   number <framenum> (0, 1, or 2). If <start_codon> is <true> a frame has to
   start with a start codon, otherwise a frame can start everywhere (i.e., at
   the first amino acid or after a stop codon). If <final_stop_codon> is <true>
   the last ORF must end with a stop codon, otherwise it can be ``open''. For
   each ORF the <orf_processor> function is called and <data>, <framenum> and
   <frame> is passed along. If <framepos> is true, the ORF range is reported in
   the coordinate system of the frame (i.e., the amino acids). Otherwise the
   coordinate system of the original sequence is used (i.e., the nucleotides).
   The correct <framenum> is needed for the conversion. */
void gt_determine_ORFs(GtORFProcessor orf_processor, void *data,
                       unsigned int framenum, const char *frame,
                       unsigned long framelen, bool start_codon,
                       bool final_stop_codon, bool framepos,
                       const char *start_codons);

#endif
