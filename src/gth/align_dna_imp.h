/*
  Copyright (c) 2003-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef ALIGN_DNA_IMP_H
#define ALIGN_DNA_IMP_H

#include "gth/align_dna.h"

#define DNA_NUMOFSCORETABLES  2

/*
  Here comes the dirty bitvector stuff to quarter the backtrace matrix.
  In this case the Retrace type defined above is only used in the lower E_STATE.
*/

#define LOWER_E_STATE_MASK      0x7     /* |0000|0111| */
#define LOWER_I_STATE_MASK      0x8     /* |0000|1000| */
#define UPPER_E_STATE_MASK      0x70    /* |0111|0000| */
#define UPPER_I_STATE_MASK      0x80    /* |1000|0000| */

#define UPPER_E_N               (DNA_E_N << 4)
#define UPPER_E_M               (DNA_E_M << 4)

#define I_STATE_E_N             (DNA_E_N << 3)
#define I_STATE_I_N             (DNA_I_N << 3)

#define UPPER_I_STATE_I_N       (I_STATE_I_N << 4)

/* these definitions refer to the state index for the two N x M score matrices
*/
typedef enum {
  DNA_E_STATE = 0,
  DNA_I_STATE,
  DNA_NUMOFSTATES
} DnaStates;

/*
  These definitions are used to retrace the optimal path in align_dna.c.
  E, I refer to the original state, N, M, and NM refer to which index
  should be decremented.

  IMPORTANT: Definition has to be consistent with dna_retracenames in
  align_dna.c.
*/
typedef enum {
  DNA_I_N = 0, /* deletion     |01|0^14
                  or intron    |01|length */
  DNA_E_N,     /* deletion     |01|0^14
                  or intron    |01|length */
  DNA_E_NM,    /* match        |00|length
                  or mismatch  |11|0^14   */
  DNA_I_NM,    /* match        |00|length
                  or mismatch  |11|0^14   */
  DNA_E_M,     /* insertion    |10|0^14   */
  DNA_I_M,     /* insertion    |10|0^14   */
  DNA_NUMOFRETRACE
} DnaRetrace;

/* the following structure bundles all tables involved in the dynamic
   programming for cDNAs/ESTs */
struct GthDPMatrix {
  GthFlt *score[DNA_NUMOFSTATES][DNA_NUMOFSCORETABLES]; /* table to store the
                                                           score of a path */
  GthPath **path;                   /* backtrace table of size
                                        gen_dp_length * ref_dp_length */
  GthPath **path_jt;
  unsigned long *intronstart[DNA_NUMOFSCORETABLES],
                *exonstart[DNA_NUMOFSCORETABLES],
                gen_dp_length,
                ref_dp_length;
};

#endif
