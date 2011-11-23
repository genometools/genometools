/*
  Copyright (c) 2004-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef ALIGN_PROTEIN_IMP_H
#define ALIGN_PROTEIN_IMP_H

#include "gth/align_protein.h"

#define WSIZE_PROTEIN   20
#define WSIZE_DNA       60 /* (3 * WSIZE_PROTEIN) */

#define GENOMICDPSTART          3
#define REFERENCEDPSTART        1

#define PROTEIN_NUMOFSCORETABLES  4

/* these definitions refer to the state index for the two N x M score matrices
*/
typedef enum {
  E_STATE = 0,
  IA_STATE,
  IB_STATE,
  IC_STATE,
  PROTEIN_NUMOFSTATES
} States;

/*
   These definitions are used to retrace the optimal path in align_protein.c.
   E, I refer to the original state, N, M, and NM refer to which index
   should be decremented.

   IMPORTANT: Definition has to be consistent with retracenames in
   align_protein.c.
*/
typedef enum {
  /*         GS2:       representation as multi editoperation: */

  E_N3M,  /* C_N3M      match               |00|00|length                   */
          /*            or mismatch         |11|00|0^12                     */
  E_N2M,  /* C_N2M      mismatch + 1 gap    |11|01|0^12                     */
  E_N1M,  /* C_N1M      mismatch + 2 gaps   |11|10|0^12                     */
  E_M,    /* C_M        insertion           |10|00|0^12                     */
  E_N3,   /* C_N3       deletion            |01|00|0^12                     */
          /*            or intron                                           */
  E_N2,   /* C_N2       deletion + 1 gap    |01|01|0^12                     */
          /*            or intron                                           */
  E_N1,   /* C_N        deletion + 2 gaps   |01|10|0^12                     */
          /*            or intron                                           */
  IA_N3M, /* I_N3M      match               |00|00|length                   */
          /*            or mismatch         |11|00|0^12                     */
  IB_N2M, /* I_N2M      match               |00|01|length (XXX: necessary?) */
          /*            or mismatch         |11|01|0^12                     */
  IC_N1M, /* I_N1M      match               |00|10|length (XXX: necessary?) */
          /*            or mismatch         |11|10|0^12                     */
  IA_N1,  /* IA_N       intron              |01|00|length                   */
  IB_N1,  /* IB_N       intron keep 1 base  |01|01|length                   */
  IC_N1,  /* IC_N       intron keep 2 bases |01|10|length                   */

  NUMOFRETRACE
} Retrace;

/* the following structure bundles core all tables involved in the DP */
typedef struct {
  /* table to store the score of a path */
  GthFlt *score[PROTEIN_NUMOFSTATES][PROTEIN_NUMOFSCORETABLES];
  GthPath **path; /* backtrace table of size gen_dp_length * ref_dp_length */
} DPtablecore;

/* structure of a path matrix byte:

   87654321
    |||+---E_STATE (bits 1 to 4)
    |||
    ||+----IA_STATE (bit 5)
    ||
    |+-----IB_STATE (bit 6)
    |
    +------IC_STATE (bit 7)

    bit 8 is unused.
*/

#define E_STATE_MASK      0xf     /* |0000|1111| */
#define IA_STATE_MASK     0x10    /* |0001|0000| */
#define IB_STATE_MASK     0x20    /* |0010|0000| */
#define IC_STATE_MASK     0x40    /* |0100|0000| */

/* the following structure bundles all tables involved in the dynamic
   programming for proteins */
struct GthDPtables {
  DPtablecore core; /* the core tables */
  unsigned long *intronstart_A[PROTEIN_NUMOFSCORETABLES],
                *intronstart_B[PROTEIN_NUMOFSCORETABLES],
                *intronstart_C[PROTEIN_NUMOFSCORETABLES];
  unsigned long *exonstart[PROTEIN_NUMOFSCORETABLES];
  unsigned char *splitcodon_B[PROTEIN_NUMOFSCORETABLES],
                *splitcodon_C1[PROTEIN_NUMOFSCORETABLES],
                *splitcodon_C2[PROTEIN_NUMOFSCORETABLES];
};

#endif
