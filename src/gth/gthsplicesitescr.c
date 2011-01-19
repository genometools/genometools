/*
  Copyright (c) 2003-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/chardef.h"
#include "core/minmax.h"
#include "core/unused_api.h"
#include "gth/gthsplicesitescr.h"

#define SPLICE_SITE_SCORE_WINDOW        50 /* (GS2 = WSIZE) */
#define SSSWINDOW_MINSIZE_FACTOR        .8

typedef struct {
  bool breaktraversealignment;
  const unsigned char *gen_seq_tran,
                      *ref_seq_tran;
  GtAlphabet *gen_alphabet;
  GthDPOptionsEST *dp_options_est;
  GthFlt splicesiteweight,    /* (GS2 = sscr) */
         maxsplicesiteweight; /* (GS2 = osscr) */
  unsigned long processedalignmentpositions;
} Calcsplicesitescoredata;

static void calcsplicesitescoreprocmismatchordeletion(Traversealignmentstate
                                                      *state,
                                                      void *data,
                                                      GT_UNUSED
                                                      unsigned long lengthofeop)
{
  Calcsplicesitescoredata *d = (Calcsplicesitescoredata*) data;
  unsigned int gen_alphabet_mapsize = gt_alphabet_size(d->gen_alphabet);
  GthDPOptionsEST *dp_options_est = d->dp_options_est;
  unsigned char genomicchar;
  if (d->processedalignmentpositions < SPLICE_SITE_SCORE_WINDOW) {
    genomicchar   = d->gen_seq_tran[state->genomicptr];
    ADDOUTPUTWEIGHTIDENTITY(d->maxsplicesiteweight, genomicchar);
    d->processedalignmentpositions++;
  }
  else
    d->breaktraversealignment = true;
}

static void calcsplicesitescoreprocinsertion(Traversealignmentstate *state,
                                             void *data,
                                             GT_UNUSED
                                             unsigned long lengthofeop)
{
  Calcsplicesitescoredata *d = (Calcsplicesitescoredata*) data;
  unsigned int gen_alphabet_mapsize = gt_alphabet_size(d->gen_alphabet);
  GthDPOptionsEST *dp_options_est = d->dp_options_est;
  unsigned char genomicchar, referencechar;
  if (d->processedalignmentpositions < SPLICE_SITE_SCORE_WINDOW) {
    genomicchar   = (unsigned char) DASH;
    referencechar = d->ref_seq_tran[state->referenceptr];
    ADDOUTPUTWEIGHT(d->splicesiteweight, genomicchar, referencechar);
    ADDOUTPUTWEIGHTIDENTITY(d->maxsplicesiteweight, genomicchar);
    d->processedalignmentpositions++;
  }
  else
    d->breaktraversealignment = true;
}

static void calcsplicesitescoreprocmatch(Traversealignmentstate *state,
                                         void *data, unsigned long lengthofeop)
{
  Calcsplicesitescoredata *d = (Calcsplicesitescoredata*) data;
  unsigned int gen_alphabet_mapsize = gt_alphabet_size(d->gen_alphabet);
  GthDPOptionsEST *dp_options_est = d->dp_options_est;
  unsigned char genomicchar, referencechar;
  unsigned long numofmatchestoprocess, alignmentpositionsleft;
  GthFlt genomicinterimvalue   = 0.0,
         referenceinterimvalue = 0.0;
  if (d->processedalignmentpositions < SPLICE_SITE_SCORE_WINDOW) {
    alignmentpositionsleft = SPLICE_SITE_SCORE_WINDOW -
                             d->processedalignmentpositions;
    numofmatchestoprocess  = MIN(lengthofeop, alignmentpositionsleft);

    genomicchar   = d->gen_seq_tran[state->genomicptr];
    referencechar = d->ref_seq_tran[state->referenceptr];
    ADDOUTPUTWEIGHT(referenceinterimvalue, genomicchar, referencechar);
    ADDOUTPUTWEIGHTIDENTITY(genomicinterimvalue, genomicchar);
    genomicinterimvalue   *= numofmatchestoprocess;
    referenceinterimvalue *= numofmatchestoprocess;
    d->splicesiteweight  += referenceinterimvalue;
    d->maxsplicesiteweight += genomicinterimvalue;
    d->processedalignmentpositions += numofmatchestoprocess;
  }
  else
   d->breaktraversealignment = true;
}

static void calcsplicesitescoreprocintron(GT_UNUSED
                                          Traversealignmentstate *state,
                                          void *data,
                                          GT_UNUSED
                                          unsigned long lengthofeop)
{
  Calcsplicesitescoredata *d = (Calcsplicesitescoredata*) data;
  d->breaktraversealignment = true;
}

static bool calcsplicesitescorebreakcondition(void *data)
{
  Calcsplicesitescoredata *d = (Calcsplicesitescoredata*) data;
  return d->breaktraversealignment;
}

void gthcalcsplicesitescore(GthDbl *splicesitescore,
                            Traversealignmentstate *oldstate,
                            const unsigned char *gen_seq_tran,
                            const unsigned char *ref_seq_tran,
                            GtAlphabet *gen_alphabet,
                            GthDPOptionsEST *dp_options_est,
                            bool acceptorsite)
{
  Traversealignmentfunctions travfunctions;
  Traversealignmentstate newstate;
  Calcsplicesitescoredata data;

  gt_assert(dp_options_est);

  travfunctions.processmismatch  = calcsplicesitescoreprocmismatchordeletion;
  travfunctions.processdeletion  = calcsplicesitescoreprocmismatchordeletion;
  travfunctions.processinsertion = calcsplicesitescoreprocinsertion;
  travfunctions.processmatch     = calcsplicesitescoreprocmatch;
  travfunctions.processintron    = calcsplicesitescoreprocintron;
  travfunctions.breakcondition   = calcsplicesitescorebreakcondition;

  /* to prevent manipulation of oldstate we copy it to newstate */
  newstate = *oldstate;

  if (!acceptorsite) { /* i.e. we want to process a donorsite */
    /* to go to the eopptr before the oldstate->eopptr */
    newstate.eopptr++;

    /* adjusting the sequence pointers to be able to correctly go backwards in
       the alignment */
    newstate.genomicptr--;
    newstate.referenceptr--;
  }

  data.breaktraversealignment      = false;
  data.gen_seq_tran                = gen_seq_tran;
  data.ref_seq_tran                = ref_seq_tran;
  data.gen_alphabet                = gen_alphabet;
  data.dp_options_est              = dp_options_est;
  data.splicesiteweight            = (GthFlt) 0.0;
  data.maxsplicesiteweight         = (GthFlt) 0.0;
  data.processedalignmentpositions = 0;

  /* for acceptorsites going forward, for donorsites going backward */
  gthtraversealignment(acceptorsite, &newstate, false, &data, &travfunctions);

  if ((data.processedalignmentpositions >=  (SSSWINDOW_MINSIZE_FACTOR *
                                                  SPLICE_SITE_SCORE_WINDOW)) &&
      (data.splicesiteweight > 0.0) &&  /* the weights must be positive */
      (data.maxsplicesiteweight > 0.0)) {
    *splicesitescore = (GthDbl) (data.splicesiteweight /
                                 data.maxsplicesiteweight);
  }
  else
    *splicesitescore = 0.0;
}
