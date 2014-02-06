/*
  Copyright (c) 2013 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

#ifndef SEQABSTRACT_H
#define SEQABSTRACT_H

#include "core/types_api.h"
#include "core/encseq_api.h"

/* Class <GtSeqabstract> represents short substrings of either <GtEncseq> or
   <GtUchar>-Arrays.
   All indices given in the methods of this class ar relative to <offset>. */
typedef struct GtSeqabstract GtSeqabstract;

GtSeqabstract* gt_seqabstract_new_empty(void);

/* Creates new <GtSeqabstract> object from <string>, starting at <offset> with
   length <len>, <string> should be long enough and <offset> within <string>.
   Ownership of <string> stays with the caller. */
GtSeqabstract* gt_seqabstract_new_gtuchar(const GtUchar *string,
                                          GtUword len,
                                          GtUword offset);

/* Creates new <GtSeqabstract> object from <encseq>, starting at <offset> with
   length <len>, fails if <offset> is out of bounds, or <offset> + <len> extends
   <encseq>.
 */
GtSeqabstract* gt_seqabstract_new_encseq(const GtEncseq *encseq,
                                         GtUword len,
                                         GtUword offset);

/* reinitialize <sa> with <string> starting at <offset> with length <len> */
void           gt_seqabstract_reinit_gtuchar(GtSeqabstract *sa,
                                             const GtUchar *string,
                                             GtUword len,
                                             GtUword offset);

/* reinitialize <sa> with <encseq> starting at <offset> with length <len> */
void           gt_seqabstract_reinit_encseq(GtSeqabstract *sa,
                                            const GtEncseq *encseq,
                                            GtUword len,
                                            GtUword offset);

/* return the length of <sa> */
GtUword        gt_seqabstract_length(const GtSeqabstract *sa);

/* return character at positon <idx> (relative to <offset>) of <sa> */
GtUchar        gt_seqabstract_encoded_char(const GtSeqabstract *sa,
                                           GtUword idx);

/* calculate longest common prefix for suffixes <ustart> and <vstart> of <useq>
   and <vseq>. */
GtUword        gt_seqabstract_lcp(bool forward,
                                  const GtSeqabstract *useq,
                                  const GtSeqabstract *vseq,
                                  GtUword ustart,
                                  GtUword vstart);

void           gt_seqabstract_delete(GtSeqabstract *sa);

#endif
