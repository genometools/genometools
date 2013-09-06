/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#ifndef SEQNUMRELPOS_H
#define SEQNUMRELPOS_H

typedef struct GtSeqnumrelpos GtSeqnumrelpos;

GtSeqnumrelpos *gt_seqnumrelpos_new(unsigned int bitsforrelpos,
                                    const GtEncseq *encseq);

void gt_seqnumrelpos_delete(GtSeqnumrelpos *snrp);

GtUword gt_seqnumrelpos_decode_pos(const GtSeqnumrelpos *snrp,
                                         GtUword seqnumrelpos);

GtUword gt_seqnumrelpos_decode_seqnum(const GtSeqnumrelpos *snrp,
                                            GtUword seqnumrelpos);

GtUword gt_seqnumrelpos_decode_relpos(const GtSeqnumrelpos *snrp,
                                            GtUword seqnumrelpos);

GtUword gt_seqnumrelpos_encode(const GtSeqnumrelpos *snrp,
                                     GtUword seqnum,
                                     GtUword relpos);

#endif
