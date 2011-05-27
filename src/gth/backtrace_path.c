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

#include "core/undef_api.h"
#include "core/safearith.h"
#include "gth/backtrace_path.h"
#include "gth/gthcutoffsminimal.h"
#include "gth/gthcutoffsrelaxed.h"
#include "gth/gthcutoffsstrict.h"

#define DETERMINE_TIMES_MAXLEN     \
        while (length >  maxlen) { \
          times_maxlen++;          \
          length -= maxlen;        \
        }

/* The following (pseudo) edit operation is used to store a dummy which is
   later replaced with a proper edit operation. */
#define DUMMY_EOP ((Editoperation) (15 << 12)) /* 11|11|0^12 */

typedef struct {
  Cutoffs start, /* the leading cutoffs */
          end;   /* the terminal cutoffs */
} Alignmentcutoffs;

struct GthBacktracePath {
  GtArray *editoperations;
                         /* assertions for an array of multi editoperations:
                            - an intron can only be of on type (00 or 01 or 11)
                            - after an 00 intron a 00 match/mismatch must follow
                            - after an 01 intron a 01 match/mismatch must follow
                            - after an 10 intron a 10 match/mismatch must follow
                         */
  GthAlphatype alphatype;
  Editoperation max_identical_length; /* defines the maximal length for a match.
                                         for dna MAXIDENTICALLENGTH should be
                                         used and for protein
                                         MAXIDENTICALLENGTH_PROTEIN */
  unsigned long gen_dp_start,
                gen_dp_length,
                ref_dp_start,
                ref_dp_length,
                dummy_index;
  Alignmentcutoffs cutoffs;
};

GthBacktracePath* gth_backtrace_path_new(unsigned long gen_dp_start,
                                         unsigned long gen_dp_length,
                                         unsigned long ref_dp_start,
                                         unsigned long ref_dp_length)
{
  GthBacktracePath *bp;
  bp = gt_calloc(1, sizeof *bp);
  bp->editoperations = gt_array_new(sizeof (Editoperation));
  bp->alphatype = UNDEF_ALPHA;
  bp->gen_dp_start = gen_dp_start;
  bp->gen_dp_length = gen_dp_length;
  bp->ref_dp_start = ref_dp_start;
  bp->ref_dp_length = ref_dp_length;
  bp->dummy_index = GT_UNDEF_ULONG;
  return bp;
}

unsigned long gth_backtrace_path_gen_dp_start(const GthBacktracePath *bp)
{
  gt_assert(bp && bp->gen_dp_start != GT_UNDEF_ULONG);
  return bp->gen_dp_start;
}

void gth_backtrace_path_set_gen_dp_start(GthBacktracePath *bp,
                                         unsigned long gen_dp_start)
{
  gt_assert(bp && gen_dp_start != GT_UNDEF_ULONG);
  bp->gen_dp_start = gen_dp_start;
}

unsigned long gth_backtrace_path_gen_dp_length(const GthBacktracePath *bp)
{
  gt_assert(bp && bp->gen_dp_length != GT_UNDEF_ULONG);
  return bp->gen_dp_length;
}

void gth_backtrace_path_set_gen_dp_length(GthBacktracePath *bp,
                                          unsigned long gen_dp_length)
{
  gt_assert(bp && gen_dp_length != GT_UNDEF_ULONG);
  bp->gen_dp_length = gen_dp_length;
}

unsigned long gth_backtrace_path_ref_dp_length(const GthBacktracePath *bp)
{
  gt_assert(bp && bp->ref_dp_length != GT_UNDEF_LONG);
  return bp->ref_dp_length;
}

void gth_backtrace_path_set_ref_dp_length(GthBacktracePath *bp,
                                          unsigned long ref_dp_length)
{
  gt_assert(bp && ref_dp_length != GT_UNDEF_ULONG);
  bp->ref_dp_length = ref_dp_length;
}

unsigned long gt_compute_indelcount(Editoperation *alignment,
                                 unsigned long alignmentlength, bool proteineop)
{
  unsigned long i, eoplength, indelcount = 0;
  Eoptype eoptype;

  for (i = 0; i < alignmentlength; i++) {
    eoptype   = gt_editoperation_type(alignment[i], proteineop);
    eoplength = gt_editoperation_length(alignment[i], proteineop);

    switch (eoptype) {
      case EOP_TYPE_MATCH:
        /* nothing to do */
        break;
      case EOP_TYPE_INTRON:
      case EOP_TYPE_INTRON_WITH_1_BASE_LEFT:
      case EOP_TYPE_INTRON_WITH_2_BASES_LEFT:
        indelcount += eoplength;
        break;
      case EOP_TYPE_MISMATCH:
        /* nothing to do */
        break;
      case EOP_TYPE_DELETION:
      case EOP_TYPE_INSERTION:
        if (proteineop)
          indelcount += eoplength * 3;
        else
          indelcount += eoplength;
        break;
      case EOP_TYPE_MISMATCH_WITH_1_GAP:
        gt_assert(proteineop);
        indelcount += eoplength;
        break;
      case EOP_TYPE_MISMATCH_WITH_2_GAPS:
        gt_assert(proteineop);
        indelcount += eoplength * 2;
        break;
      case EOP_TYPE_DELETION_WITH_1_GAP:
      case EOP_TYPE_DELETION_WITH_2_GAPS:
        if (proteineop)
          indelcount += eoplength * 3;
        else
          indelcount += eoplength;
        break;
      default: gt_assert(0);
    }
  }

  return indelcount;
}

unsigned long gth_backtrace_path_indelcount(const GthBacktracePath *bp)
{
  gt_assert(bp);
  gt_assert(bp->alphatype == DNA_ALPHA || bp->alphatype == PROTEIN_ALPHA);
  return gt_compute_indelcount(gt_array_get_space(bp->editoperations),
                            gt_array_size(bp->editoperations),
                            bp->alphatype == PROTEIN_ALPHA);
}

unsigned long gth_backtrace_path_genomiccutoff_start(const GthBacktracePath *bt)
{
  gt_assert(bt);
  return bt->cutoffs.start.genomiccutoff;
}

unsigned long gth_backtrace_path_referencecutoff_start(const
                                                       GthBacktracePath *bt)
{
  gt_assert(bt);
  return bt->cutoffs.start.referencecutoff;
}

unsigned long gth_backtrace_path_eopcutoff_start(const GthBacktracePath *bt)
{
  gt_assert(bt);
  return bt->cutoffs.start.eopcutoff;
}

unsigned long gth_backtrace_path_genomiccutoff_end(const GthBacktracePath *bt)
{
  gt_assert(bt);
  return bt->cutoffs.end.genomiccutoff;
}

unsigned long gth_backtrace_path_referencecutoff_end(const GthBacktracePath *bt)
{
  gt_assert(bt);
  return bt->cutoffs.end.referencecutoff;
}

unsigned long gth_backtrace_path_eopcutoff_end(const GthBacktracePath *bt)
{
  gt_assert(bt);
  return bt->cutoffs.end.eopcutoff;
}

void gth_backtrace_path_set_cutoffs_start(GthBacktracePath *bt,
                                          Cutoffs *cutoffs)
{
  gt_assert(bt && cutoffs);
  bt->cutoffs.start = *cutoffs;
}

void gth_backtrace_path_set_cutoffs_end(GthBacktracePath *bt, Cutoffs *cutoffs)
{
  gt_assert(bt && cutoffs);
  bt->cutoffs.end = *cutoffs;
}

GthAlphatype gth_backtrace_path_alphatype(const GthBacktracePath *bp)
{
  gt_assert(bp);
  gt_assert(bp->alphatype == DNA_ALPHA || bp->alphatype == PROTEIN_ALPHA);
  return bp->alphatype;
}

void gth_backtrace_path_set_alphatype(GthBacktracePath *bp,
                                      GthAlphatype alphatype)
{
  gt_assert(bp && bp->alphatype == UNDEF_ALPHA);
  gt_assert(alphatype == DNA_ALPHA || alphatype == PROTEIN_ALPHA);
  bp->alphatype = alphatype;
  if (alphatype == DNA_ALPHA)
    bp->max_identical_length = MAXIDENTICALLENGTH;
  else
    bp->max_identical_length = MAXIDENTICALLENGTH_PROTEIN;
}

static void determine_cutoffs(GthBacktracePath *bp,
                              GthCutoffmode leadcutoffsmode,
                              GthCutoffmode termcutoffsmode,
                              unsigned long cutoffsminexonlen)
{
  Traversealignmentfunctions travfunctions;
  Traversealignmentstate travstate;
  Relaxedcutoffsdata relaxedcutoffsdata;
  Strictcutoffsdata strictcutoffsdata;
  Minimalcutoffsdata minimalcutoffsdata;
  bool proteineop = bp->alphatype == PROTEIN_ALPHA;

  /* sum of edit operations equals referencelength (before cutoffs)", */
  gt_assert(gth_backtrace_path_is_valid(bp));

  /* setting the traverse alignment state */
  travstate.proteineop      = proteineop;
  travstate.processing_intron_with_1_base_left  = false;
  travstate.processing_intron_with_2_bases_left = false;
  travstate.alignment       = gth_backtrace_path_get(bp);
  travstate.alignmentlength = gth_backtrace_path_length(bp);
  travstate.eopptr          = travstate.alignment +
                              travstate.alignmentlength - 1;
  travstate.genomicptr      = 0;
  travstate.referenceptr    = 0;

  /* cutting of leading indels in the sequences */
  switch (leadcutoffsmode) {
    case RELAXED:
      gt_initRelaxedcutoffsTravfunctions(&travfunctions);
      gt_initRelaxedcutoffsdata(&relaxedcutoffsdata, &bp->cutoffs.start);
      gthtraversealignment(true, &travstate, proteineop, &relaxedcutoffsdata,
                           &travfunctions);
      break;
    case STRICT:
      gt_initStrictcutoffsTravfunctions(&travfunctions);
      gt_initStrictcutoffsdata(&strictcutoffsdata, &bp->cutoffs.start,
                            cutoffsminexonlen);
      gthtraversealignment(true , &travstate , proteineop, &strictcutoffsdata,
                           &travfunctions);
      break;
    case MINIMAL:
      gt_initMinimalcutoffsTravfunctions(&travfunctions);
      gt_initMinimalcutoffsdata(&minimalcutoffsdata, &bp->cutoffs.start);
      gthtraversealignment(true, &travstate, proteineop, &minimalcutoffsdata,
                           &travfunctions);
      break;
    default: gt_assert(0);
  }

  /* resetting the traverse alignment state */
  travstate.processing_intron_with_1_base_left  = false;
  travstate.processing_intron_with_2_bases_left = false;
  travstate.eopptr = gth_backtrace_path_get(bp);
  travstate.genomicptr = 0;
  travstate.referenceptr = 0;

  /* cutting of terminal indels in the sequences */
  switch (termcutoffsmode) {
    case RELAXED:
      gt_initRelaxedcutoffsTravfunctions(&travfunctions);
      gt_initRelaxedcutoffsdata(&relaxedcutoffsdata, &bp->cutoffs.end);
      gthtraversealignment(false, &travstate, proteineop, &relaxedcutoffsdata,
                           &travfunctions);
      break;
    case STRICT:
      gt_initStrictcutoffsTravfunctions(&travfunctions);
      gt_initStrictcutoffsdata(&strictcutoffsdata, &bp->cutoffs.end,
                            cutoffsminexonlen);
      gthtraversealignment(false, &travstate, proteineop, &strictcutoffsdata,
                           &travfunctions);
      break;
    case MINIMAL:
      gt_initMinimalcutoffsTravfunctions(&travfunctions);
      gt_initMinimalcutoffsdata(&minimalcutoffsdata, &bp->cutoffs.end);
      gthtraversealignment(false, &travstate, proteineop, &minimalcutoffsdata,
                           &travfunctions);
      break;
    default: gt_assert(0);
  }

  /* sum of edit operations equals referencelength (after cutoffs) */
  gt_assert(gth_backtrace_path_is_valid(bp));
}

void gth_backtrace_path_determine_cutoffs(GthBacktracePath *bp,
                                          GthCutoffmode leadcutoffsmode,
                                          GthCutoffmode termcutoffsmode,
                                          unsigned long cutoffsminexonlen)
{
  gt_assert(bp);
  /* ref_alphatype is valid */
  gt_assert(bp->alphatype == DNA_ALPHA || bp->alphatype == PROTEIN_ALPHA);
  memset(&bp->cutoffs, 0, sizeof bp->cutoffs);
  determine_cutoffs(bp, leadcutoffsmode, termcutoffsmode, cutoffsminexonlen);
}

static bool is_insertion(Editoperation eop)
{
  if (eop == INSERTIONEOP)
    return true;
  return false;
}

static bool is_intron(Editoperation eop)
{
  if ((eop & ~MAXIDENTICALLENGTH) == DELETIONEOP &&
      (eop & MAXIDENTICALLENGTH) > 0) {
    return true;
  }
  return false;
}

/* The following function removes zero base exons in <alignment>.
   That is, a continuous stretch of insertions between two introns is moved past
   the intron(s) to the right.
   This is necessary for the gthcomputescores() function to work correctly.
   Otherwise one would get exons with borders (i, i-1) which can not be
   processed by the successive functions. */
void gt_remove_zero_base_exons(Editoperation *alignment, long alignmentlength,
                            GthStat *stat)
{
  long i, j;
  for (i = 1; i < alignmentlength - 1; i++) {
    if (is_insertion(alignment[i])) {
      /* editoperation is insertion, check for surrounding introns */
      if (is_intron(alignment[i-1])) {
        /* intron to the left -> go to the right */
        for (j = i + 1; j < alignmentlength; j++) {
          if (is_intron(alignment[j])) {
            /* Insertion(s) surrounded by Introns -> move insertion(s) past the
               (complete) intron */
            while (j < alignmentlength) {
              Editoperation tmp_eop;

              /* swap */
              tmp_eop = alignment[i];
              alignment[i] = alignment[j];
              alignment[j] = tmp_eop;
              i++;
              j++;
              if (!is_intron(alignment[j]))
                break;
            }
            /* increase counter */
            gth_stat_increment_numofremovedzerobaseexons(stat);
            break;
          }
          else if (!is_insertion(alignment[j]))
            break; /* Insertions are not surrounded by Introns */
        }
      }
    }
  }
}

void gth_backtrace_path_remove_zero_base_exons(GthBacktracePath *bp,
                                               GthStat *stat)
{
  gt_assert(bp);
  gt_remove_zero_base_exons(gth_backtrace_path_get(bp),
                         gth_backtrace_path_length(bp), stat);
}

static bool contains_no_zero_base_exons(Editoperation *alignment,
                                        long alignmentlength)
{
  long i, j;

  for (i = 1; i < alignmentlength - 1; i++) {
    if (is_insertion(alignment[i])) {
      /* editoperation is Insertion, check for surrounding Introns */
      if (is_intron(alignment[i-1])) {
        /* intron to the left -> go to the right */
        for (j = i + 1; j < alignmentlength; j++) {
          if (is_intron(alignment[j])) {
            /* Insertion is surrounded by Introns -> return false */
            return false;
          }          else if (!is_insertion(alignment[j]))
            break; /* Insertions are not surrounded by Introns */
        }      }
    }
  }
  /* no zero base exons -> return true */
  return true;
}

bool gth_backtrace_path_contains_no_zero_base_exons(const GthBacktracePath *bp)
{
  gt_assert(bp);
  return contains_no_zero_base_exons(gth_backtrace_path_get(bp) +
                                     bp->cutoffs.end.eopcutoff,
                                     gth_backtrace_path_length(bp) -
                                     bp->cutoffs.start.eopcutoff -
                                     bp->cutoffs.end.eopcutoff);
}

static void add_eop_type_to_eop_array(GtArray *bp, Eoptype eoptype,
                                      unsigned long length, bool proteineop)
{
  Editoperation eop,
                maxlen = proteineop ? (Editoperation) MAXIDENTICALLENGTH_PROTEIN
                                    : (Editoperation) MAXIDENTICALLENGTH;
  Eoptype tmp_eoptype;
  unsigned long i, times_maxlen = 0;

  gt_assert(length > 0);

  switch (eoptype) {
    case EOP_TYPE_MATCH:
      /* here we reproduce the artifact resulting from the dummys used in the
         backtracing procedure to make sure that the parsed array of edit
         operations is exactly the same as the one we have in memory */
      if (proteineop && /* this needs only to be checked for protein bp */
          length > 1 &&       /* and when the length is larger 1 */
          gt_array_size(bp)) { /* we have already stored an eop */
        tmp_eoptype = gt_editoperation_type(*(Editoperation*)
                                         gt_array_get_last(bp), proteineop);
        if (tmp_eoptype == EOP_TYPE_INTRON_WITH_1_BASE_LEFT ||
            tmp_eoptype == EOP_TYPE_INTRON_WITH_2_BASES_LEFT) {
          eop = 1;
          gt_array_add(bp, eop);
          length--;
        }
      }

      /* we store the eop which has not maximal length first to make sure that
         after reversing the array of editoperations has the same form as the
         original one */
      DETERMINE_TIMES_MAXLEN;
      gt_assert(length > 0);
      eop = (Editoperation) length;
      gt_array_add(bp, eop);
      for (i = 0; i < times_maxlen; i++)
        gt_array_add(bp, maxlen);
      break;
    case EOP_TYPE_INTRON:
      DETERMINE_TIMES_MAXLEN;
      eop  = DELETIONEOP;
      eop += length;
      gt_array_add(bp, eop);
      eop  = DELETIONEOP;
      eop += maxlen;
      for (i = 0; i < times_maxlen; i++)
        gt_array_add(bp, eop);
      break;
    case EOP_TYPE_INTRON_WITH_1_BASE_LEFT:
      DETERMINE_TIMES_MAXLEN;
      eop  = DELETION_WITH_1_GAP_EOP;
      eop += length;
      gt_array_add(bp, eop);
      eop  = DELETION_WITH_1_GAP_EOP;
      eop += maxlen;
      for (i = 0; i < times_maxlen; i++)
        gt_array_add(bp, eop);
      break;
    case EOP_TYPE_INTRON_WITH_2_BASES_LEFT:
      DETERMINE_TIMES_MAXLEN;
      eop  = DELETION_WITH_2_GAPS_EOP;
      eop += length;
      gt_array_add(bp, eop);
      eop  = DELETION_WITH_2_GAPS_EOP;
      eop += maxlen;
      for (i = 0; i < times_maxlen; i++)
        gt_array_add(bp, eop);
      break;
    case EOP_TYPE_MISMATCH:
      eop = MISMATCHEOP;
      for (i = 0; i < length; i++)
        gt_array_add(bp, eop);
      break;
    case EOP_TYPE_DELETION:
      eop = DELETIONEOP;
      for (i = 0; i < length; i++)
        gt_array_add(bp, eop);
      break;
    case EOP_TYPE_INSERTION:
      eop = INSERTIONEOP;
      for (i = 0; i < length; i++)
        gt_array_add(bp, eop);
      break;
    case EOP_TYPE_MISMATCH_WITH_1_GAP:
      eop = MISMATCH_WITH_1_GAP_EOP;
      for (i = 0; i < length; i++)
        gt_array_add(bp, eop);
      break;
    case EOP_TYPE_MISMATCH_WITH_2_GAPS:
      eop = MISMATCH_WITH_2_GAPS_EOP;
      for (i = 0; i < length; i++)
        gt_array_add(bp, eop);
      break;
    case EOP_TYPE_DELETION_WITH_1_GAP:
      eop = DELETION_WITH_1_GAP_EOP;
      for (i = 0; i < length; i++)
        gt_array_add(bp, eop);
      break;
    case EOP_TYPE_DELETION_WITH_2_GAPS:
      eop = DELETION_WITH_2_GAPS_EOP;
      for (i = 0; i < length; i++)
        gt_array_add(bp, eop);
      break;
    default: gt_assert(0);
  }
}

void gth_backtrace_path_add_eop(GthBacktracePath *bp, Eoptype eoptype,
                            unsigned long length)
{
  gt_assert(bp);
  gt_assert(bp->alphatype == DNA_ALPHA || bp->alphatype == PROTEIN_ALPHA);
  add_eop_type_to_eop_array(bp->editoperations, eoptype, length,
                            bp->alphatype == PROTEIN_ALPHA);
}

void gth_backtrace_path_add_match(GthBacktracePath *bp,
                                  bool ensure_single_match)
{
  Editoperation *eopptr, match_eop = 1;
  unsigned long eopid, lenid;
  gt_assert(bp);
  gt_assert(bp->alphatype == DNA_ALPHA || bp->alphatype == PROTEIN_ALPHA);
  if (!gt_array_size(bp->editoperations) || ensure_single_match)
    gt_array_add(bp->editoperations, match_eop);
  else {
    eopptr = gt_array_get_last(bp->editoperations);
    eopid  = *eopptr & ~bp->max_identical_length;
    lenid  = *eopptr &  bp->max_identical_length;
    if (eopid == 0 && lenid > 0 && lenid < bp->max_identical_length)
      (*eopptr)++;
    else
      gt_array_add(bp->editoperations, match_eop);
  }
}

void gth_backtrace_path_add_mismatch(GthBacktracePath *bp)
{
  Editoperation mismatch_eop = MISMATCHEOP;
  gt_assert(bp);
  gt_assert(bp->alphatype == DNA_ALPHA || bp->alphatype == PROTEIN_ALPHA);
  gt_array_add(bp->editoperations, mismatch_eop);
}

void gth_backtrace_path_add_deletion(GthBacktracePath *bp)
{
  Editoperation deletion_eop = DELETIONEOP;
  gt_assert(bp);
  gt_assert(bp->alphatype == DNA_ALPHA || bp->alphatype == PROTEIN_ALPHA);
  gt_array_add(bp->editoperations, deletion_eop);
}

void gth_backtrace_path_add_insertion(GthBacktracePath *bp)
{
  Editoperation insertion_eop = INSERTIONEOP;
  gt_assert(bp);
  gt_assert(bp->alphatype == DNA_ALPHA || bp->alphatype == PROTEIN_ALPHA);
  gt_array_add(bp->editoperations, insertion_eop);
}

void gth_backtrace_path_add_intron(GthBacktracePath *bp)
{
  Editoperation *eopptr, intron_eop = DELETIONEOP + 1;
  unsigned long eopid, lenid;
  gt_assert(bp);
  gt_assert(bp->alphatype == DNA_ALPHA || bp->alphatype == PROTEIN_ALPHA);
  if (!gt_array_size(bp->editoperations))
    gt_array_add(bp->editoperations, intron_eop);
  else {
    eopptr = gt_array_get_last(bp->editoperations);
    eopid  = *eopptr & ~bp->max_identical_length;
    lenid  = *eopptr &  bp->max_identical_length;
    if (eopid == DELETIONEOP && lenid > 0 && lenid < bp->max_identical_length)
      (*eopptr)++;
    else
      gt_array_add(bp->editoperations, intron_eop);
  }
}

void gth_backtrace_path_add_mismatch_with_1_gap(GthBacktracePath *bp)
{
  Editoperation mismatch_with_1_gap_eop = MISMATCH_WITH_1_GAP_EOP;
  gt_assert(bp);
  gt_assert(bp->alphatype == PROTEIN_ALPHA);
  gt_assert(bp->max_identical_length == MAXIDENTICALLENGTH_PROTEIN);
  gt_array_add(bp->editoperations, mismatch_with_1_gap_eop);
}

void gth_backtrace_path_add_mismatch_with_2_gaps(GthBacktracePath *bp)
{
  Editoperation mismatch_with_2_gaps_eop = MISMATCH_WITH_2_GAPS_EOP;
  gt_assert(bp);
  gt_assert(bp->alphatype == PROTEIN_ALPHA);
  gt_assert(bp->max_identical_length == MAXIDENTICALLENGTH_PROTEIN);
  gt_array_add(bp->editoperations, mismatch_with_2_gaps_eop);
}

void gth_backtrace_path_add_deletion_with_1_gap(GthBacktracePath *bp)
{
  Editoperation deletion_with_1_gap_eop = DELETION_WITH_1_GAP_EOP;
  gt_assert(bp);
  gt_assert(bp->alphatype == PROTEIN_ALPHA);
  gt_assert(bp->max_identical_length == MAXIDENTICALLENGTH_PROTEIN);
  gt_array_add(bp->editoperations, deletion_with_1_gap_eop);
}

void gth_backtrace_path_add_deletion_with_2_gaps(GthBacktracePath *bp)
{
  Editoperation deletion_with_2_gaps_eop = DELETION_WITH_2_GAPS_EOP;
  gt_assert(bp);
  gt_assert(bp->alphatype == PROTEIN_ALPHA);
  gt_assert(bp->max_identical_length == MAXIDENTICALLENGTH_PROTEIN);
  gt_array_add(bp->editoperations, deletion_with_2_gaps_eop);
}

void gth_backtrace_path_add_intron_with_1_base_left(GthBacktracePath *bp)
{
  Editoperation *eopptr,
                intron_with_1_base_left_eop = DELETION_WITH_1_GAP_EOP + 1;
  unsigned long eopid, lenid;
  gt_assert(bp);
  gt_assert(bp->alphatype == PROTEIN_ALPHA);
  gt_assert(bp->max_identical_length == MAXIDENTICALLENGTH_PROTEIN);
  if (!gt_array_size(bp->editoperations))
    gt_array_add(bp->editoperations, intron_with_1_base_left_eop);
  else {
    eopptr = gt_array_get_last(bp->editoperations);
    eopid  = *eopptr & ~bp->max_identical_length;
    lenid  = *eopptr &  bp->max_identical_length;
    if (eopid ==  DELETION_WITH_1_GAP_EOP && lenid > 0 &&
        lenid < bp->max_identical_length) {
      (*eopptr)++;
    }
    else
      gt_array_add(bp->editoperations, intron_with_1_base_left_eop);
  }
}

void gth_backtrace_path_add_intron_with_2_bases_left(GthBacktracePath *bp)
{
  Editoperation *eopptr,
                intron_with_2_bases_left_eop = DELETION_WITH_2_GAPS_EOP + 1;
  unsigned long eopid, lenid;
  gt_assert(bp);
  gt_assert(bp->alphatype == PROTEIN_ALPHA);
  gt_assert(bp->max_identical_length == MAXIDENTICALLENGTH_PROTEIN);
  if (!gt_array_size(bp->editoperations))
    gt_array_add(bp->editoperations, intron_with_2_bases_left_eop);
  else {
    eopptr = gt_array_get_last(bp->editoperations);
    eopid  = *eopptr & ~bp->max_identical_length;
    lenid  = *eopptr &  bp->max_identical_length;
    if (eopid ==  DELETION_WITH_2_GAPS_EOP && lenid > 0 &&
        lenid < bp->max_identical_length) {
      (*eopptr)++;
    }
    else
      gt_array_add(bp->editoperations, intron_with_2_bases_left_eop);
  }
}

void gth_backtrace_path_add_dummy(GthBacktracePath *bp)
{
  Editoperation dummy_eop = DUMMY_EOP;
  gt_assert(bp);
  gt_assert(bp->alphatype == PROTEIN_ALPHA);
  gt_assert(bp->max_identical_length == MAXIDENTICALLENGTH_PROTEIN);
  gt_assert(bp->dummy_index == GT_UNDEF_ULONG);
  gt_array_add(bp->editoperations, dummy_eop);
  bp->dummy_index = gt_array_size(bp->editoperations) - 1;
}

void gth_backtrace_path_set_dummy(GthBacktracePath *bp, bool match)
{
  Editoperation *eopptr;
  gt_assert(bp);
  gt_assert(bp->alphatype == PROTEIN_ALPHA);
  gt_assert(bp->max_identical_length == MAXIDENTICALLENGTH_PROTEIN);
  gt_assert(bp->dummy_index != GT_UNDEF_ULONG);
  eopptr = gt_array_get(bp->editoperations, bp->dummy_index);
  if (match)
    *eopptr = 1; /* match */
  else
    *eopptr = MISMATCHEOP;
  bp->dummy_index = GT_UNDEF_ULONG; /* reset dummy */
}

bool gth_backtrace_path_contain_dummy(const GthBacktracePath *bp)
{
  gt_assert(bp);
  if (bp->dummy_index == GT_UNDEF_ULONG)
    return false;
  return true;
}

bool gth_backtrace_path_last_is_intron(const GthBacktracePath *bp)
{
  Eoptype eoptype;

  gt_assert(bp);

  /* check if a dummy has just been inserted */
  if (bp->dummy_index != GT_UNDEF_ULONG &&
      gt_array_size(bp->editoperations) - 1 == bp->dummy_index) {
    return false;
  }

  eoptype = gt_editoperation_type(*(Editoperation*)
                               gt_array_get_last(bp->editoperations),
                               bp->alphatype == PROTEIN_ALPHA);
  if (eoptype == EOP_TYPE_INTRON ||
      eoptype == EOP_TYPE_INTRON_WITH_1_BASE_LEFT ||
      eoptype == EOP_TYPE_INTRON_WITH_2_BASES_LEFT) {
    return true;
  }
  return false;

}

void gth_backtrace_path_reverse(GthBacktracePath *bp)
{
  Editoperation *front, *back, tmp;
  unsigned long i = 0;
  gt_assert(bp);
  for (front = gt_array_get_space(bp->editoperations),
       back  = (Editoperation*) gt_array_get_space(bp->editoperations) +
               gt_array_size(bp->editoperations) - 1;
       front < back; front++, back--, i++) {
    tmp = *front;
    *front = *back;
    *back = tmp;
  }
}

static void ensure_eop_of_len_1_before_introns(GtArray *editoperations)
{
  Editoperation eop, *eopptr;
  Eoptype eoptype;
  unsigned long eoplength;
  GtArray *backup;
  bool processing_necessary = false,
       split_match          = false;

  /* check if processing is necessary
     the check is rather simple, it might be possible that
     ``processing_necessary'' is set to ``true'' whereas in fact no processing
     is necessary */
  for (eopptr = gt_array_get_space(editoperations);
       eopptr < (Editoperation*) gt_array_get_space(editoperations) +
                                 gt_array_size(editoperations) - 1;
       eopptr++) {
    if ((eoptype = gt_editoperation_type(*eopptr, true)) ==
        EOP_TYPE_INTRON_WITH_1_BASE_LEFT ||
        eoptype == EOP_TYPE_INTRON_WITH_2_BASES_LEFT) {
      processing_necessary = true;
      break;
    }
  }

  if (processing_necessary) {
    /* init backup for the editoperations */
    backup = gt_array_new(sizeof (Editoperation));

    /* fill backup */
    gt_array_add_array(backup, editoperations);

    /* reset the original edit operations */
    gt_array_set_size(editoperations, 0);

    /* process the backup and fill the original editoperations */
    for (eopptr = gt_array_get_space(backup);
         eopptr < (Editoperation*)
                  gt_array_get_space(backup) + gt_array_size(backup);
         eopptr++) {

      if ((eoptype = gt_editoperation_length(*eopptr, true)) ==
          EOP_TYPE_INTRON_WITH_1_BASE_LEFT ||
          eoptype == EOP_TYPE_INTRON_WITH_2_BASES_LEFT) {
        split_match = true;
      }
      else if (split_match) {
        if (eoptype == EOP_TYPE_MATCH) {
          split_match = false;
          if ((eoplength = gt_editoperation_length(*eopptr, true)) > 1) {
            eop = 1;
            gt_array_add(editoperations, eop);
            eop = eoplength - 1;
            gt_array_add(editoperations, eop);
            continue;
          }
        }
        else if (eoptype == EOP_TYPE_MISMATCH ||
                 eoptype == EOP_TYPE_MISMATCH_WITH_1_GAP) {
          split_match = false;
        }
      }
      gt_array_add(editoperations, *eopptr);
    }

    /* free backup */
    gt_array_delete(backup);
  }
}

void gth_backtrace_path_ensure_length_1_before_introns(GthBacktracePath *bp)
{
  if (bp->alphatype == PROTEIN_ALPHA)
    ensure_eop_of_len_1_before_introns(bp->editoperations);
}

size_t gth_backtrace_path_sizeof (const GthBacktracePath *bp)
{
  gt_assert(bp);
  /* XXX -> array_sizeof () */
  return gt_array_size(bp->editoperations) * sizeof (Editoperation) +
         sizeof (GthBacktracePath);
}

void gth_backtrace_path_show(const GthBacktracePath *bp, bool xmlout,
                             unsigned int indentlevel, GtFile *outfp)
{
  gt_assert(bp);
  gt_assert(bp->alphatype == DNA_ALPHA || bp->alphatype == PROTEIN_ALPHA);
  gt_editoperation_show(gth_backtrace_path_get(bp),
                        gth_backtrace_path_length(bp),
                        bp->alphatype == PROTEIN_ALPHA, xmlout, indentlevel,
                        outfp);
}

void gth_backtrace_path_show_complete(const GthBacktracePath *bp, bool xmlout,
                                      unsigned int indentlevel,
                                      GtFile *outfp)
{
  gt_assert(bp);
  gt_assert(bp->alphatype == DNA_ALPHA || bp->alphatype == PROTEIN_ALPHA);
  gt_editoperation_show(gt_array_get_space(bp->editoperations),
                     gt_array_size(bp->editoperations),
                     bp->alphatype == PROTEIN_ALPHA, xmlout, indentlevel,
                     outfp);
}

void gth_backtrace_path_cutoff_start(GthBacktracePath *bp)
{
  gt_assert(bp);
  gt_assert(bp->gen_dp_start != GT_UNDEF_ULONG);
  gt_assert(bp->gen_dp_length != GT_UNDEF_ULONG);
  gt_assert(bp->ref_dp_start != GT_UNDEF_ULONG);
  gt_assert(bp->ref_dp_length != GT_UNDEF_ULONG);
  if (bp->cutoffs.start.genomiccutoff) {
    bp->gen_dp_start += bp->cutoffs.start.genomiccutoff;
    gt_assert(bp->gen_dp_length >= bp->cutoffs.start.genomiccutoff);
    bp->gen_dp_length -= bp->cutoffs.start.genomiccutoff;
    bp->cutoffs.start.genomiccutoff = 0;
  }
  if (bp->cutoffs.start.referencecutoff) {
    bp->ref_dp_start += bp->cutoffs.start.referencecutoff;
    gt_assert(bp->ref_dp_length >= bp->cutoffs.start.referencecutoff);
    bp->ref_dp_length -= bp->cutoffs.start.referencecutoff;
    bp->cutoffs.start.referencecutoff = 0;
  }
  if (bp->cutoffs.start.eopcutoff) {
    gt_array_set_size(bp->editoperations,
                      gt_array_size(bp->editoperations) -
                      bp->cutoffs.start.eopcutoff);
    bp->cutoffs.start.eopcutoff = 0;
  }
}

void gth_backtrace_path_cutoff_end(GthBacktracePath *bp)
{
  gt_assert(bp);
  gt_assert(bp->gen_dp_length != GT_UNDEF_ULONG);
  gt_assert(bp->ref_dp_length != GT_UNDEF_ULONG);
  if (bp->cutoffs.end.genomiccutoff) {
    gt_assert(bp->gen_dp_length >= bp->cutoffs.end.genomiccutoff);
    bp->gen_dp_length -= bp->cutoffs.end.genomiccutoff;
    bp->cutoffs.end.genomiccutoff = 0;
  }
  if (bp->cutoffs.end.referencecutoff) {
    gt_assert(bp->ref_dp_length >= bp->cutoffs.end.referencecutoff);
    bp->ref_dp_length -= bp->cutoffs.end.referencecutoff;
    bp->cutoffs.end.referencecutoff = 0;
  }
  if (bp->cutoffs.end.eopcutoff) {
    gt_array_rem_span(bp->editoperations, 0, bp->cutoffs.end.eopcutoff - 1);
    bp->cutoffs.end.eopcutoff = 0;
  }
}

#ifndef NDEBUG
static bool backtrace_path_start_cutoffs_are_set(const GthBacktracePath *bp)
{
  gt_assert(bp);
  if (bp->cutoffs.start.genomiccutoff   ||
      bp->cutoffs.start.referencecutoff ||
      bp->cutoffs.start.eopcutoff) {
    return true;
  }
  return false;
}

static bool backtrace_path_end_cutoffs_are_set(const GthBacktracePath *bp)
{
  gt_assert(bp);
  if (bp->cutoffs.end.genomiccutoff   ||
      bp->cutoffs.end.referencecutoff ||
      bp->cutoffs.end.eopcutoff) {
    return true;
  }
  return false;
}

static bool backtrace_path_cutoffs_are_set(const GthBacktracePath *bp)
{
  gt_assert(bp);
  if (backtrace_path_start_cutoffs_are_set(bp) ||
      backtrace_path_end_cutoffs_are_set(bp)) {
    return true;
  }
  return false;
}
#endif

void gth_backtrace_path_cutoff_walked_path(GthBacktracePath *bp,
                                           const GthPathWalker *pw,
                                           bool showeops, GtFile *outfp)
{
  unsigned int length;
  gt_assert(bp && pw);
  if (gth_path_walker_is_forward(pw)) {
    gt_assert(!backtrace_path_start_cutoffs_are_set(bp));
    if (showeops) {
      gt_file_xprintf(outfp, "%s(): show path walker\n", __func__);
      gth_path_walker_show(pw, outfp);
      gt_file_xprintf(outfp, "%s(): show backtrace path (before eop "
                         "removal)\n", __func__);
      gth_backtrace_path_show(bp, false, 0, outfp);
    }
    /* remove complete eops */
    gt_array_set_size(bp->editoperations,
                      gt_array_size(bp->editoperations) -
                      gth_path_walker_actual_eops(pw));
    if (showeops) {
      gt_file_xprintf(outfp, "%s(): show backtrace path (after eop "
                         "removal)\n", __func__);
      gth_backtrace_path_show(bp, false, 0, outfp);
    }
    /* remove part of last eop */
    if (gth_path_walker_steps_in_current_eop(pw)) {
      length = gt_editoperation_length(*(Editoperation*)
                                    gt_array_get_last(bp->editoperations),
                                    bp->alphatype == PROTEIN_ALPHA);
      gt_assert(length > gth_path_walker_steps_in_current_eop(pw));
      gt_editoperation_set_length(gt_array_get_last(bp->editoperations),
                               length-gth_path_walker_steps_in_current_eop(pw),
                               bp->alphatype == PROTEIN_ALPHA);
    }
    /* adjusting genomic and reference DP ranges */
    bp->gen_dp_start += gth_path_walker_gen_distance(pw);
    bp->gen_dp_length -= gth_path_walker_gen_distance(pw);
    bp->ref_dp_start += gth_path_walker_ref_distance(pw);
    bp->ref_dp_length -= gth_path_walker_ref_distance(pw);
  }
  else {
    gt_assert(0); /* XXX: implement reverse case */
    gt_assert(!backtrace_path_end_cutoffs_are_set(bp));
  }
}

/* XXX: remove */
#if 0
static void cutoff_end_refseq(GthBacktracePath *bp, unsigned long reflength)
{

  unsigned long eoplength, i = 0;
  bool breakloop = false;
  Editoperation *eop;
  Eoptype eoptype;
  gt_assert(bp && reflength);
  gt_assert(bp->alphatype == DNA_ALPHA || bp->alphatype == PROTEIN_ALPHA);

  for (;;) {
    eop = (Editoperation*) gt_array_get(bp->editoperations, i);
    eoptype   = gt_editoperation_type(*eop, bp->alphatype == PROTEIN_ALPHA);
    eoplength = gt_editoperation_length(*eop, bp->alphatype == PROTEIN_ALPHA);
    i++;

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
        if (eoplength >= reflength) {
          breakloop = true;
          if (eoplength > reflength) {
            gt_assert(eoplength > 2);
            *eop &= ~bp->max_identical_length;
            *eop |= eoplength - 1;
            i--;
          }
        }
        break;
      default: gt_assert(0);
    }
    if (breakloop)
      break;
    reflength -= eoplength;
  }

  if (i)
    gt_array_rem_span(bp->editoperations, 0, i-1);
}
#endif

void gth_backtrace_path_prepend(GthBacktracePath *out_bp,
                                const GthBacktracePath *in_bp)
{
  gt_assert(out_bp && in_bp);
  gt_assert(!backtrace_path_start_cutoffs_are_set(out_bp));
  gt_assert(!backtrace_path_cutoffs_are_set(in_bp));
  gt_assert(in_bp->gen_dp_start + in_bp->gen_dp_length == out_bp->gen_dp_start);
  gt_assert(in_bp->ref_dp_start + in_bp->ref_dp_length == out_bp->ref_dp_start);
  out_bp->gen_dp_start = in_bp->gen_dp_start;
  out_bp->gen_dp_length += in_bp->gen_dp_length;
  out_bp->ref_dp_start = in_bp->ref_dp_start;
  out_bp->ref_dp_length += in_bp->ref_dp_length;
  gt_array_add_array(out_bp->editoperations, in_bp->editoperations);
}

void gth_backtrace_path_append(GthBacktracePath *out_bp,
                               const GthBacktracePath *in_bp)
{
  gt_assert(out_bp && in_bp);
  gt_assert(!backtrace_path_end_cutoffs_are_set(out_bp));
  gt_assert(!backtrace_path_cutoffs_are_set(in_bp));
  gt_assert(out_bp->gen_dp_start+out_bp->gen_dp_length == in_bp->gen_dp_start);
  gt_assert(out_bp->ref_dp_start+out_bp->ref_dp_length == in_bp->ref_dp_start);
  out_bp->gen_dp_length += in_bp->gen_dp_length;
  out_bp->ref_dp_length += in_bp->ref_dp_length;
  gt_array_prepend_array(out_bp->editoperations, in_bp->editoperations);
}

bool gth_backtrace_path_is_valid(const GthBacktracePath *bp)
{
  bool is_valid;
  gt_assert(bp);
  gt_assert(bp->alphatype == DNA_ALPHA || bp->alphatype == PROTEIN_ALPHA);
  gt_assert(bp->ref_dp_length != GT_UNDEF_ULONG);
  is_valid =
    gt_eops_equal_referencelength((Editoperation*)
                               gt_array_get_space(bp->editoperations)
                               + bp->cutoffs.end.eopcutoff,
                               gt_safe_cast2long(gt_array_size(bp
                                                              ->editoperations))
                               - bp->cutoffs.start.eopcutoff
                               - bp->cutoffs.end.eopcutoff,
                               gt_safe_cast2long(bp->ref_dp_length)
                               - bp->cutoffs.start.referencecutoff
                               - bp->cutoffs.end.referencecutoff,
                               bp->alphatype == PROTEIN_ALPHA);
  return is_valid;
}

void gth_backtrace_path_reset(GthBacktracePath *bp)
{
  gt_assert(bp);
  gt_array_set_size(bp->editoperations, 0);
  bp->alphatype = UNDEF_ALPHA;
  bp->dummy_index = GT_UNDEF_ULONG;
}

void gth_backtrace_path_delete(GthBacktracePath *bp)
{
  if (!bp) return;
  gt_array_delete(bp->editoperations);
  gt_free(bp);
}

Editoperation* gth_backtrace_path_get(const GthBacktracePath *bp)
{
  gt_assert(bp);
  return (Editoperation*) gt_array_get_space(bp->editoperations)
         + bp->cutoffs.end.eopcutoff;
}

unsigned long gth_backtrace_path_length(const GthBacktracePath *bp)
{
  unsigned long num_of_eops;
  gt_assert(bp);
  num_of_eops = gt_array_size(bp->editoperations);
  if (bp->cutoffs.start.eopcutoff + bp->cutoffs.end.eopcutoff >= num_of_eops)
    return 0; /* prevent overflow */
  return num_of_eops - bp->cutoffs.start.eopcutoff - bp->cutoffs.end.eopcutoff;
}
