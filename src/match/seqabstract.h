/*
  Copyright (c) 2015 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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
   All indices given in the methods of this class ar relative to <startpos>. */
typedef struct GtSeqabstract GtSeqabstract;

GtSeqabstract* gt_seqabstract_new_empty(void);

/* Reset <GtSeqabstract> object to initial values */
void gt_seqabstract_reset(GtSeqabstract *sa);

/* Creates new <GtSeqabstract> object from <string>, starting at <startpos> with
   length <len>, <string> should be long enough and <startpos> within <string>.
   Ownership of <string> stays with the caller. */
GtSeqabstract* gt_seqabstract_new_gtuchar(bool rightextension,
                                          GtReadmode readmode,
                                          const GtUchar *string,
                                          GtUword len,
                                          GtUword startpos,
                                          GtUword totallength);

/* Creates new <GtSeqabstract> object from <encseq>, starting at <startpos> with
   length <len>, fails if <startpos> is out of bounds, or
   <startpos> + <len> extends <encseq>.
 */
GtSeqabstract* gt_seqabstract_new_encseq(bool rightextension,
                                         GtReadmode readmode,
                                         const GtEncseq *encseq,
                                         GtUword len,
                                         GtUword startpos);

/* reinitialize <sa> with <string> starting at <startpos> with length <len> */
void           gt_seqabstract_reinit_gtuchar(bool rightextension,
                                             GtReadmode readmode,
                                             GtSeqabstract *sa,
                                             const GtUchar *string,
                                             GtUword len,
                                             GtUword startpos,
                                             GtUword totallength);

/* reinitialize <sa> with <encseq> starting at <startpos> with length <len> */
void           gt_seqabstract_reinit_encseq(bool rightextension,
                                            GtReadmode readmode,
                                            GtSeqabstract *sa,
                                            const GtEncseq *encseq,
                                            GtUword len,
                                            GtUword startpos);

/* return the length of <sa> */
GtUword        gt_seqabstract_length(const GtSeqabstract *sa);

/* return character at positon <idx> (relative to <startpos>) of <sa> */
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

/* set the readmode flag which is GT_READMODE_FORWARD by default */

void gt_seqabstract_readmode_set(GtSeqabstract *sa,GtReadmode readmode);

char *gt_seqabstract_get(bool rightextension,const GtSeqabstract *seq);

void gt_seqabstract_seqstartpos_set(GtSeqabstract *sa,GtUword seqstartpos);

void gt_seqabstract_totallength_set(GtSeqabstract *sa,GtUword totallength);

#endif
