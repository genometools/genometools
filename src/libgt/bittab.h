/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef BITTAB_H
#define BITTAB_H

#include <stdbool.h>
#include <stdio.h>
#include "array.h"

typedef struct Bittab Bittab;

Bittab*       bittab_new(unsigned long num_of_bits);
void          bittab_set_bit(Bittab*, unsigned long);
void          bittab_unset_bit(Bittab*, unsigned long);
void          bittab_complement(Bittab*, const Bittab*);           /* a=~b   */
void          bittab_equal(Bittab*, const Bittab*);                /* a=b    */
void          bittab_and(Bittab*, const Bittab*, const Bittab*);   /* a=b&c  */
void          bittab_or(Bittab*, const Bittab*, const Bittab*);    /* a=b|c  */
void          bittab_nand(Bittab*, const Bittab*, const Bittab*);  /* a=b&~c */
void          bittab_and_equal(Bittab*, const Bittab*);            /* a&=b   */
void          bittab_or_equal(Bittab*, const Bittab*);             /* a|=b   */
void          bittab_shift_left_equal(Bittab*);                    /* a<<=1  */
void          bittab_shift_right_equal(Bittab*);                   /* a>>=1  */
void          bittab_unset(Bittab*);                               /* a=0    */
void          bittab_show(const Bittab*, FILE*);
void          bittab_get_all_bitnums(const Bittab*, Array*);
bool          bittab_bit_is_set(const Bittab*, unsigned long);
bool          bittab_is_true(const Bittab*);
bool          bittab_cmp(const Bittab*, const Bittab*);
unsigned long bittab_get_first_bitnum(const Bittab*);
unsigned long bittab_get_last_bitnum(const Bittab*);
unsigned long bittab_get_next_bitnum(const Bittab*, unsigned long);
unsigned long bittab_count_set_bits(Bittab*);
unsigned long bittab_size(Bittab*);
int           bittab_unit_test(Error*);
void          bittab_delete(Bittab*);

#if 0
  a typical iterator loop:

  for (bit  = bittab_get_first_bitnum(b);
       bit != bittab_get_last_bitnum(b);
       bit  = bittab_get_next_bitnum(b, bit)) {
    /* ... */
  }
#endif

#endif
