/*
  Copyright (c) 2003-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef GTHTRAVALIGN_H
#define GTHTRAVALIGN_H

#include <stdbool.h>
#include "gth/gthalignment.h"

typedef struct {
  bool proteineop,             /* processing protein edit operations */
       processing_intron_with_1_base_left,
       processing_intron_with_2_bases_left;
  Editoperation *alignment,    /* the editoperations of the whole alignment and
                                  (doesn't change when using
                                   gthtraversealignment()) */
                *eopptr;       /* the pointer to the current editoperation of
                                  the alignment (changes) */
  long alignmentlength,        /* its length, i.e. the number of Editoperations
                                */
       genomicptr,             /* the pointer to the current character of the
                                  genomic sequence */
       referenceptr,           /* the pointer to the current character of the
                                  reference sequence */
       firstbaseleftptr,
       secondbaseleftptr;
} Traversealignmentstate;

typedef struct {
  void (*processmismatch)(Traversealignmentstate*, void*, unsigned long);
  void (*processdeletion)(Traversealignmentstate*, void*, unsigned long);
  void (*processinsertion)(Traversealignmentstate*, void*, unsigned long);
  void (*processmatch)(Traversealignmentstate*, void*, unsigned long);
  void (*processintron)(Traversealignmentstate*, void*, unsigned long);
  bool(*breakcondition)(void*);
  /* additional functions for protein edit operations */
  void (*processintron_with_1_base_left)(Traversealignmentstate*, void*,
                                         unsigned long);
  void (*processintron_with_2_bases_left)(Traversealignmentstate*, void*,
                                          unsigned long);
  void (*processmismatch_with_1_gap)(Traversealignmentstate*, void*,
                                     unsigned long);
  void (*processmismatch_with_2_gaps)(Traversealignmentstate*, void*,
                                      unsigned long);
  void (*processdeletion_with_1_gap)(Traversealignmentstate*, void*,
                                     unsigned long);
  void (*processdeletion_with_2_gaps)(Traversealignmentstate*, void*,
                                      unsigned long);
} Traversealignmentfunctions;

/* The following function checks if the sum of multi edit operations
   (sum in terms of reference sequence bases) equals the reference length.
   If alignmentlength is less than 0 it is true returned, no matter what
   referencelength is. */
bool gt_eops_equal_referencelength(Editoperation *alignment,
                                   long alignmentlength,
                                   long referencelength, bool proteineop);

/* The following function traverses an alignment given by an
   Traversealignmentstate structure <state> in forward direction, if
   <forward> equals true, in backward direction otherwise.
   Forward direction means forward in the alignment (i.e. from left to right),
   but backward (i.e. from right to left) in the Array of Editoperations.
   For every Editoperation visited during the traversion, the corresponding
   Traversealignmentfunction is called, if it is defined.
   Thereby, <proteineop> denotes if protein edit operations are used. */
void gthtraversealignment(bool forward, Traversealignmentstate *state,
                          bool proteineop, void *data,
                          Traversealignmentfunctions *travfunctions);

#endif
