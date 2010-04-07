/*
  Copyright (c) 2003-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2003-2007 Center for Bioinformatics, University of Hamburg

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

#ifndef EDITOPERATION_H
#define EDITOPERATION_H

#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include "core/file.h"

/*
  An edit operation is an uint16_t value which stores
  in the two most significant bits one of the following edit operation
  or 0 in which case the remaining 14 bits are used for
  representing the length of a pair of identical substrings in the alignment.
  The largest such pair of identical substrings to store is thus of length
  MAXIDENTICALLENGTH.

  An intron is represented as a DELETIONEOP where the remaining 14 bits
  represent the length (01|length).
*/

typedef uint16_t Editoperation;

#define MAXIDENTICALLENGTH ((1 << 14) - 1)
#define DELETIONEOP        ((Editoperation) (1 << 14)) /* 01|00|0^12 */
#define INSERTIONEOP       ((Editoperation) (1 << 15)) /* 10|00|0^12 */
#define MISMATCHEOP        ((Editoperation) (3 << 14)) /* 11|00|0^12 */

/*
  To be able to represent DNA/protein alignments we need additional multi edit
  operations. We use two additional bits for this purpose (the third and
  fourth most significant ones).
  In this case the largest pair of identical substrings or introns which can
  be stored in a multi edit operation shrinks to MAXIDENTICALLENGTH_PROTEIN.

  The introns used in DNA alignments (01|00|length) represent introns
  which start after a completely processed codon in DNA/protein alignments.
  Introns which start after an incompletely processed codon are represented as
  follows:
  - Introns which start after two bases of a codon are represented as a
    DELETION_WITH_1_GAP_EOP, because 1 base is left.
    The remaining 12 bits represent the length (01|01|length).
  - Introns wich start after one base of a codon are represented as a
    DELETION_WITH_2_GAPS_EOP, because 2 bases are left.
    The remaining 12 bits represent the length (01|10|length).
*/

#define MAXIDENTICALLENGTH_PROTEIN ((1 << 12) - 1)
#define MISMATCH_WITH_1_GAP_EOP  ((Editoperation) (13 << 12)) /* 11|01|0^12 */
#define MISMATCH_WITH_2_GAPS_EOP ((Editoperation) (14 << 12)) /* 11|10|0^12 */
#define DELETION_WITH_1_GAP_EOP  ((Editoperation) (5 << 12))  /* 01|01|0^12 */
#define DELETION_WITH_2_GAPS_EOP ((Editoperation) (6 << 12))  /* 01|10|0^12 */

typedef enum {
  EOP_TYPE_MATCH = 0,
  EOP_TYPE_INTRON,
  EOP_TYPE_INTRON_WITH_1_BASE_LEFT,
  EOP_TYPE_INTRON_WITH_2_BASES_LEFT,
  EOP_TYPE_MISMATCH,
  EOP_TYPE_DELETION,
  EOP_TYPE_INSERTION,
  EOP_TYPE_MISMATCH_WITH_1_GAP,
  EOP_TYPE_MISMATCH_WITH_2_GAPS,
  EOP_TYPE_DELETION_WITH_1_GAP,
  EOP_TYPE_DELETION_WITH_2_GAPS,
  NUM_OF_EOP_TYPES
} Eoptype;

Eoptype      gt_editoperation_type(Editoperation, bool proteineop);
unsigned int gt_editoperation_length(Editoperation, bool proteineop);
void         gt_editoperation_set_length(Editoperation*, unsigned int,
                                         bool proteineop);

void         gt_editoperation_show(Editoperation *eops,
                                   unsigned long num_of_eops,
                                   bool proteineops, bool xmlout,
                                   unsigned int indentlevel,
                                   GtFile *outfp);

#endif
