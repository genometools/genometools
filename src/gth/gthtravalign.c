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

#include "core/codon_api.h"
#include "gth/gthtravalign.h"

bool gt_eops_equal_referencelength(Editoperation *alignment,
                                long alignmentlength,
                                long referencelength,
                                bool proteineop)
{
  Eoptype eoptype;
  unsigned long eoplength;
  long i, sumofeops = 0;

  for (i = 0; i < alignmentlength; i++) {
    eoptype   = gt_editoperation_type(alignment[i], proteineop);
    eoplength = gt_editoperation_length(alignment[i], proteineop);

    switch (eoptype) {
      case EOP_TYPE_DELETION:
      case EOP_TYPE_DELETION_WITH_1_GAP:
      case EOP_TYPE_DELETION_WITH_2_GAPS:
      case EOP_TYPE_INTRON:
      case EOP_TYPE_INTRON_WITH_1_BASE_LEFT:
      case EOP_TYPE_INTRON_WITH_2_BASES_LEFT:
        /* nothing to do */
        break;
      case EOP_TYPE_INSERTION:
      case EOP_TYPE_MISMATCH:
      case EOP_TYPE_MISMATCH_WITH_1_GAP:
      case EOP_TYPE_MISMATCH_WITH_2_GAPS:
      case EOP_TYPE_MATCH:
        sumofeops += eoplength;
        break;
      default: gt_assert(0);
    }
  }

  if (alignmentlength < 0 || sumofeops == referencelength)
    return true;
  return false;
}

static void match_mismatch_after_travfunc(bool forward,
                                          Traversealignmentstate *state,
                                          bool proteineop,
                                          unsigned long eoplength)
{
  if (proteineop) {
    if (state->processing_intron_with_1_base_left) {
      /* we have to assure this for correct functioning of the STRICT cutoffs */
      gt_assert(eoplength == 1);
      state->processing_intron_with_1_base_left = false;
      if (forward) {
        state->genomicptr   += eoplength * GT_CODON_LENGTH - 1;
        state->referenceptr += eoplength;
      }
      else {
        state->genomicptr   -= eoplength * GT_CODON_LENGTH - 1;
        state->referenceptr -= eoplength;
      }
    }
    else if (state->processing_intron_with_2_bases_left) {
      /* we have to assure this for correct functioning of the STRICT cutoffs */
      gt_assert(eoplength == 1);
      state->processing_intron_with_2_bases_left = false;
      if (forward) {
        state->genomicptr   += eoplength * GT_CODON_LENGTH - 2;
        state->referenceptr += eoplength - 1;
      }
      else {
        state->genomicptr   -= eoplength * GT_CODON_LENGTH - 2;
        state->referenceptr -= eoplength - 1;
      }
    }
    else {
      if (forward) {
        state->genomicptr   += eoplength * GT_CODON_LENGTH;
        state->referenceptr += eoplength;
      }
      else {
        state->genomicptr   -= eoplength * GT_CODON_LENGTH;
        state->referenceptr -= eoplength;
      }
    }
  }
  else {
    if (forward) {
      state->genomicptr   += eoplength;
      state->referenceptr += eoplength;
    }
    else {
      state->genomicptr   -= eoplength;
      state->referenceptr -= eoplength;
    }
  }
}

void gthtraversealignment(bool forward, Traversealignmentstate *state,
                          bool proteineop, void *data,
                          Traversealignmentfunctions *travfunctions)
{
  Eoptype eoptype;
  unsigned long eoplength;

  for (;forward ? state->eopptr >= state->alignment
                : state->eopptr <= state->alignment
                                 + state->alignmentlength - 1;
        forward ? state->eopptr-- : state->eopptr++) {
    eoptype   = gt_editoperation_type(*state->eopptr, proteineop);
    eoplength = gt_editoperation_length(*state->eopptr, proteineop);

    /* we are not processing two intron types at the same time */
    gt_assert(!(state->processing_intron_with_1_base_left &&
                state->processing_intron_with_2_bases_left));

    switch (eoptype) {
      case EOP_TYPE_MISMATCH:
        if (travfunctions->processmismatch)
          travfunctions->processmismatch(state, data, eoplength);
        match_mismatch_after_travfunc(forward, state, proteineop, eoplength);
        break;
      case EOP_TYPE_DELETION:
        gt_assert(eoplength == 1);
        if (travfunctions->processdeletion)
          travfunctions->processdeletion(state, data, eoplength);
        if (proteineop) {
          if (forward)
            state->genomicptr += GT_CODON_LENGTH;
          else
            state->genomicptr -= GT_CODON_LENGTH;
        }
        else {
          if (forward)
            state->genomicptr++;
          else
            state->genomicptr--;
        }
        break;
      case EOP_TYPE_INSERTION:
        gt_assert(eoplength == 1);
        if (travfunctions->processinsertion)
          travfunctions->processinsertion(state, data, eoplength);
        if (forward)
          state->referenceptr++;
        else
          state->referenceptr--;
        break;
      case EOP_TYPE_MATCH:
        if (travfunctions->processmatch)
          travfunctions->processmatch(state, data, eoplength);
        match_mismatch_after_travfunc(forward, state, proteineop, eoplength);
        break;
      case EOP_TYPE_INTRON:
        /* we are not processing with 1 base left here */
        gt_assert(!state->processing_intron_with_1_base_left);
        /* we are not processing with 2 bases left here */
        gt_assert(!state->processing_intron_with_2_bases_left);
        if (travfunctions->processintron)
          travfunctions->processintron(state, data, eoplength);
        if (forward)
          state->genomicptr += eoplength;
        else
          state->genomicptr -= eoplength;
        break;
      case EOP_TYPE_INTRON_WITH_1_BASE_LEFT:
        gt_assert(proteineop);
        /* we are not processing with 2 bases left here */
        gt_assert(!state->processing_intron_with_2_bases_left);
        if (!state->processing_intron_with_1_base_left) {
          state->processing_intron_with_1_base_left = true;
          state->firstbaseleftptr = state->genomicptr;
          if (forward)
            state->genomicptr++;
          else
            state->genomicptr--;
        }
        if (travfunctions->processintron_with_1_base_left)
          travfunctions->processintron_with_1_base_left(state, data, eoplength);
        if (forward)
          state->genomicptr += eoplength;
        else
          state->genomicptr -= eoplength;
        break;
      case EOP_TYPE_INTRON_WITH_2_BASES_LEFT:
        gt_assert(proteineop);
        /* we are not processing with 1 base left here */
        gt_assert(!state->processing_intron_with_1_base_left);
        if (!state->processing_intron_with_2_bases_left) {
          state->processing_intron_with_2_bases_left = true;
          state->firstbaseleftptr = state->genomicptr;
          if (forward)
            state->genomicptr++;
          else
            state->genomicptr--;
          state->secondbaseleftptr = state->genomicptr;
          if (forward) {
            state->genomicptr++;
            state->referenceptr++;
          }
          else {
            state->genomicptr--;
            state->referenceptr--;
          }
        }
        if (travfunctions->processintron_with_2_bases_left) {
          travfunctions->processintron_with_2_bases_left(state, data,
                                                         eoplength);
        }
        if (forward)
          state->genomicptr += eoplength;
        else
          state->genomicptr -= eoplength;
        break;
      case EOP_TYPE_MISMATCH_WITH_1_GAP:
        gt_assert(proteineop && eoplength == 1);
        if (travfunctions->processmismatch_with_1_gap)
          travfunctions->processmismatch_with_1_gap(state, data, eoplength);
        if (state->processing_intron_with_1_base_left) {
          state->processing_intron_with_1_base_left = false;
          if (forward) {
            state->genomicptr++;
            state->referenceptr++;
          }
          else {
            state->genomicptr--;
            state->referenceptr--;
          }
        }
        else if (state->processing_intron_with_2_bases_left)
          state->processing_intron_with_2_bases_left = false;
        else {
          if (forward) {
            state->genomicptr  += 2;
            state->referenceptr++;
          }
          else {
            state->genomicptr  -= 2;
            state->referenceptr--;
          }
        }
        break;
      case EOP_TYPE_MISMATCH_WITH_2_GAPS:
        gt_assert(proteineop && eoplength == 1);
        if (travfunctions->processmismatch_with_2_gaps)
          travfunctions->processmismatch_with_2_gaps(state, data, eoplength);
        if (state->processing_intron_with_1_base_left) {
          state->processing_intron_with_1_base_left = false;
          /* the genomicptr needs no further processing in this case */
          if (forward)
            state->referenceptr++;
          else
            state->referenceptr--;
        }
        else {
          if (forward) {
            state->genomicptr++;
            state->referenceptr++;
          }
          else {
            state->genomicptr--;
            state->referenceptr--;
          }
        }
        break;
      case EOP_TYPE_DELETION_WITH_1_GAP:
        gt_assert(proteineop && eoplength == 1);
        if (travfunctions->processdeletion_with_1_gap)
          travfunctions->processdeletion_with_1_gap(state, data, eoplength);
        if (state->processing_intron_with_1_base_left) {
          state->processing_intron_with_1_base_left = false;
          if (forward)
            state->genomicptr += 1;
          else
            state->genomicptr -= 1;
        }
        else {
          if (forward)
            state->genomicptr += 2;
          else
            state->genomicptr -= 2;
        }
        break;
      case EOP_TYPE_DELETION_WITH_2_GAPS:
        gt_assert(proteineop && eoplength == 1);
        if (travfunctions->processdeletion_with_2_gaps)
          travfunctions->processdeletion_with_2_gaps(state, data, eoplength);
        if (!state->processing_intron_with_1_base_left) {
          if (forward)
            state->genomicptr++;
          else
            state->genomicptr--;
        }
        break;
      default: gt_assert(0);
    }

    if (travfunctions->breakcondition && travfunctions->breakcondition(data))
      break;
  }
}
