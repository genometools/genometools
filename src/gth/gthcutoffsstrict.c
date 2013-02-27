/*
  Copyright (c) 2003-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/codon_api.h"
#include "core/safearith.h"
#include "core/unused_api.h"
#include "gth/gthcutoffsstrict.h"

static void strictcutoffsprocindelorintron(GT_UNUSED
                                           Traversealignmentstate *state,
                                           void *data,
                                           GT_UNUSED unsigned long lengthofeop)
{
  Strictcutoffsdata *d = (Strictcutoffsdata*) data;
  d->cutoffs->eopcutoff++;
  d->actualexonlengthgenomic   = 0;
  d->actualexonlengthreference = 0;
  d->actualexonnumofeops       = 0;
}

static void strictcutoffsprocmatchormismatch_generic(Traversealignmentstate
                                                     *state, void *data,
                                                     unsigned long lengthofeop,
                                                     unsigned long
                                                     genomic_eoplen_proteins)
{
  Strictcutoffsdata *d = (Strictcutoffsdata*) data;
  unsigned long eop_length_factor = state->proteineop
                                    ? genomic_eoplen_proteins
                                    : 1;

  if (!state->processing_intron_with_1_base_left &&
      !state->processing_intron_with_2_bases_left &&
      d->actualexonlengthgenomic + lengthofeop * eop_length_factor >=
      d->cutoffsminexonlen) {
    gt_safe_sub(d->cutoffs->genomiccutoff,
                (unsigned long) gt_safe_labs(state->genomicptr),
                d->actualexonlengthgenomic);
    gt_safe_sub(d->cutoffs->referencecutoff,
                (unsigned long) gt_safe_labs(state->referenceptr),
                d->actualexonlengthreference);
    d->cutoffs->eopcutoff = d->cutoffs->eopcutoff - d->actualexonnumofeops;
    d->breakforloop = true;
  }
  else {
    d->cutoffs->eopcutoff++;
    d->actualexonlengthgenomic   += lengthofeop * eop_length_factor;
    d->actualexonlengthreference += lengthofeop;
    d->actualexonnumofeops++;
  }
}

static void strictcutoffsprocmatchormismatch(Traversealignmentstate *state,
                                             void *data,
                                             unsigned long lengthofeop)
{
  strictcutoffsprocmatchormismatch_generic(state, data, lengthofeop,
                                           GT_CODON_LENGTH);
}

static void strictcutoffsprocmismatch_with_1_gap(Traversealignmentstate *state,
                                                 void *data,
                                                 unsigned long lengthofeop)
{
  return strictcutoffsprocmatchormismatch_generic(state, data, lengthofeop, 2);
}

static void strictcutoffsprocmismatch_with_2_gaps(Traversealignmentstate *state,
                                                  void *data,
                                                  unsigned long lengthofeop)
{
  return strictcutoffsprocmatchormismatch_generic(state, data, lengthofeop, 1);
}

static bool strictcutoffsbreakcondition(void *data)
{
  Strictcutoffsdata *d = (Strictcutoffsdata*) data;
  return d->breakforloop;
}

void gt_initStrictcutoffsTravfunctions(Traversealignmentfunctions
                                       *travfunctions)
{
  travfunctions->processmismatch  = strictcutoffsprocmatchormismatch;
  travfunctions->processdeletion  = strictcutoffsprocindelorintron;
  travfunctions->processinsertion = strictcutoffsprocindelorintron;
  travfunctions->processmatch     = strictcutoffsprocmatchormismatch;
  travfunctions->processintron    = strictcutoffsprocindelorintron;
  travfunctions->breakcondition   = strictcutoffsbreakcondition;

  /* additional functions for protein edit operations */
  travfunctions->processintron_with_1_base_left  =
    strictcutoffsprocindelorintron;
  travfunctions->processintron_with_2_bases_left =
    strictcutoffsprocindelorintron;
  travfunctions->processmismatch_with_1_gap      =
    strictcutoffsprocmismatch_with_1_gap;
  travfunctions->processmismatch_with_2_gaps     =
    strictcutoffsprocmismatch_with_2_gaps;
  travfunctions->processdeletion_with_1_gap      =
    strictcutoffsprocindelorintron;
  travfunctions->processdeletion_with_2_gaps     =
    strictcutoffsprocindelorintron;
}

void gt_initStrictcutoffsdata(Strictcutoffsdata *strictcutoffsdata,
                              Cutoffs *cutoffs, unsigned long cutoffsminexonlen)
{
  strictcutoffsdata->cutoffs                   = cutoffs;
  strictcutoffsdata->breakforloop              = false;
  strictcutoffsdata->cutoffsminexonlen         = cutoffsminexonlen;
  strictcutoffsdata->actualexonlengthgenomic   = 0;
  strictcutoffsdata->actualexonlengthreference = 0;
  strictcutoffsdata->actualexonnumofeops       = 0;
}
