/*
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#ifndef SEQINFO_API_H
#define SEQINFO_API_H

/* Holds information about one sequence in a <GtEncseq>.
   The field <seqstartpos> contains the position of the first character
   in the encoded sequence while <seqlength> contains the length of the
   sequence. */
typedef struct GtSeqinfo GtSeqinfo;

struct GtSeqinfo
{
  unsigned long seqstartpos;
  unsigned long seqlength;
};

/* Creates a new <GtSeqinfo> object with both values initialized to zero. */
GtSeqinfo*    gt_seqinfo_new(void);
/* Creates a new <GtSeqinfo> object with both values initialized to <startpos>
   and <length>. */
GtSeqinfo*    gt_seqinfo_new_with_values(unsigned long startpos,
                                         unsigned long length);
/* Returns the start position stored in <si>. */
unsigned long gt_seqinfo_startpos(GtSeqinfo *si);
/* Returns the sequence length stored in <si>. */
unsigned long gt_seqinfo_length(GtSeqinfo *si);
/* Sets the start position in <si> to <startpos>. */
void          gt_seqinfo_set_startpos(GtSeqinfo *si, unsigned long startpos);
/* Sets the sequence length in <si> to <length>. */
void          gt_seqinfo_set_length(GtSeqinfo *si, unsigned long length);
/* Deletes <si>. */
void          gt_seqinfo_delete(GtSeqinfo *si);

#endif
