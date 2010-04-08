/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007-2010 Center for Bioinformatics, University of Hamburg

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

#ifndef ENCODEDSEQUENCE_INLINED_H
#define ENCODEDSEQUENCE_INLINED_H

#include "core/encodedsequence_rep.h"

#define gt_encodedsequence_total_length(ENCSEQ) \
          ((ENCSEQ)->totallength)

#define gt_encodedsequence_num_of_sequences(ENCSEQ) \
          ((ENCSEQ)->numofdbsequences)

#define GT_REVERSEPOS(TOTALLENGTH,POS) \
          ((TOTALLENGTH) - 1 - (POS))

#define GT_MAKECOMPL(CC) \
          (ISSPECIAL(CC) ? (CC) : (GtUchar) 3 - (CC))

static inline GtUchar gt_encodedsequence_get_encoded_char(
                                              const GtEncodedsequence *encseq,
                                              unsigned long pos,
                                              GtReadmode readmode)
{
  return (readmode == GT_READMODE_FORWARD)
          ? encseq->plainseq[pos]
          : ((readmode == GT_READMODE_REVERSE)
            ? encseq->plainseq[GT_REVERSEPOS(encseq->totallength,pos)]
            : ((readmode == GT_READMODE_COMPL)
              ? GT_MAKECOMPL(encseq->plainseq[pos])
              : GT_MAKECOMPL(encseq->plainseq[
                           GT_REVERSEPOS(encseq->totallength,pos)])
              )
            )
         ;
}

#define gt_encodedsequence_extract_encoded_char(ENCSEQ,POS,RM) \
          gt_encodedsequence_get_encoded_char(ENCSEQ,POS,RM)

#define gt_encodedsequence_get_encoded_charnospecial(ENCSEQ,POS,RM) \
          gt_encodedsequence_get_encoded_char(ENCSEQ,POS,RM)

#define gt_encodedsequence_get_encoded_char_sequential(ENCSEQ, \
                                                   ENCSEQSTATE,POS,READMODE) \
          gt_encodedsequence_get_encoded_char(ENCSEQ,POS,READMODE)
                                             GtEncodedsequenceScanstate *esr);

#endif
