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

#include "core/safearith.h"
#include "core/unused_api.h"
#include "gth/gthcutoffsrelaxed.h"

static void relaxedcutoffsprocmatchormismatch(Traversealignmentstate *state,
                                              void *data,
                                              GT_UNUSED
                                              unsigned long lengthofeop)
{
  Relaxedcutoffsdata *d = (Relaxedcutoffsdata*) data;

  if (state->processing_intron_with_1_base_left ||
      state->processing_intron_with_2_bases_left) {
    /* in this case we process a single match or a single mismatch after an
       intron with 1 or 2 bases left.
       we skip this single (mis)match because otherwise the pointer would
       become asynchronous and this is the best solution to solve this problem.
     */
    gt_assert(lengthofeop == 1);
    d->cutoffs->eopcutoff++;
  }
  else {
    d->cutoffs->genomiccutoff = gt_safe_labs(state->genomicptr);
    d->cutoffs->referencecutoff = gt_safe_labs(state->referenceptr);
    d->breakforloop = true;
  }
}

static void relaxedcutoffsprocindelorintron(GT_UNUSED
                                            Traversealignmentstate *state,
                                            void *data,
                                            GT_UNUSED
                                            unsigned long lengthofeop)
{
  Relaxedcutoffsdata *d = (Relaxedcutoffsdata*) data;
  d->cutoffs->eopcutoff++;
}

static bool relaxedcutoffsbreakcondition(void *data)
{
  Relaxedcutoffsdata *d = (Relaxedcutoffsdata*) data;
  return d->breakforloop;
}

void gt_initRelaxedcutoffsTravfunctions(Traversealignmentfunctions
                                        *travfunctions)
{
  travfunctions->processmismatch  = relaxedcutoffsprocmatchormismatch;
  travfunctions->processdeletion  = relaxedcutoffsprocindelorintron;
  travfunctions->processinsertion = relaxedcutoffsprocindelorintron;
  travfunctions->processmatch     = relaxedcutoffsprocmatchormismatch;
  travfunctions->processintron    = relaxedcutoffsprocindelorintron;
  travfunctions->breakcondition   = relaxedcutoffsbreakcondition;

  /* additional functions for protein edit operations */
  travfunctions->processintron_with_1_base_left  =
    relaxedcutoffsprocindelorintron;
  travfunctions->processintron_with_2_bases_left =
    relaxedcutoffsprocindelorintron;
  travfunctions->processmismatch_with_1_gap      =
    relaxedcutoffsprocmatchormismatch;
  travfunctions->processmismatch_with_2_gaps     =
    relaxedcutoffsprocmatchormismatch;
  travfunctions->processdeletion_with_1_gap      =
    relaxedcutoffsprocindelorintron;
  travfunctions->processdeletion_with_2_gaps     =
    relaxedcutoffsprocindelorintron;
}

void gt_initRelaxedcutoffsdata(Relaxedcutoffsdata *relaxedcutoffsdata,
                               Cutoffs *cutoffs)
{
  relaxedcutoffsdata->cutoffs      = cutoffs;
  relaxedcutoffsdata->breakforloop = false;
}
