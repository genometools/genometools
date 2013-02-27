/*
  Copyright (c) 2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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
#include "gth/path_walker.h"

struct GthPathWalker {
  bool forward,
       proteineop,
       processing_intron_with_1_base_left,
       processing_intron_with_2_bases_left;
  Editoperation *alignment,
                *eopptr;
  long alignmentlength;
  Eoptype last_eop_type;
  unsigned long last_eop_length,
                eop_distance,
                gen_distance,
                ref_distance,
                actual_eops;
  unsigned int  steps_in_current_eop;
};

GthPathWalker* gth_path_walker_new(const GthBacktracePath *bp, bool forward)
{
  GthPathWalker *pw;
  gt_assert(bp);
  gt_assert(forward); /* XXX: implement reverse walking */
  pw = gt_calloc(1, sizeof *pw);
  pw->forward = forward;
  pw->proteineop = gth_backtrace_path_alphatype(bp) == PROTEIN_ALPHA;
  pw->alignment = gth_backtrace_path_get(bp);
  pw->alignmentlength = gth_backtrace_path_length(bp);
  pw->eopptr = pw->alignment + pw->alignmentlength - 1;
  pw->last_eop_type = NUM_OF_EOP_TYPES;
  return pw;
}

void gth_path_walker_delete(GthPathWalker *pw)
{
  if (!pw) return;
  gt_free(pw);
}

bool gth_path_walker_is_forward(const GthPathWalker *pw)
{
  gt_assert(pw);
  return pw->forward;
}

bool gth_path_walker_has_next(const GthPathWalker *pw)
{
  gt_assert(pw);
  if ((pw->last_eop_length) ||
      (pw->forward && pw->eopptr >= pw->alignment) ||
      (!pw->forward && pw->eopptr <= pw->alignment + pw->alignmentlength - 1)) {
    return true;
  }
  return false;
}

static void step(GthPathWalker *pw)
{
  gt_assert(pw && pw->last_eop_length);

  /* we are not processing two intron types at the same time */
  gt_assert(!(pw->processing_intron_with_1_base_left &&
              pw->processing_intron_with_2_bases_left));

  pw->eop_distance++;

  switch (pw->last_eop_type) {
    case EOP_TYPE_MATCH:
    case EOP_TYPE_MISMATCH:
      if (pw->proteineop) {
        if (pw->processing_intron_with_1_base_left) {
          pw->processing_intron_with_1_base_left = false;
          pw->gen_distance += GT_CODON_LENGTH - 1;
          pw->ref_distance += 1;
        }
        else if (pw->processing_intron_with_2_bases_left) {
          pw->processing_intron_with_2_bases_left = false;
          pw->gen_distance += GT_CODON_LENGTH - 2;
        }
        else {
          pw->gen_distance += GT_CODON_LENGTH;
          pw->ref_distance += 1;
        }
      }
      else {
        pw->gen_distance++;
        pw->ref_distance++;
      }
      break;
    case EOP_TYPE_DELETION:
      if (pw->proteineop)
        pw->gen_distance += GT_CODON_LENGTH;
      else
        pw->gen_distance++;
      break;
    case EOP_TYPE_INSERTION:
      pw->ref_distance++;
      break;
    case EOP_TYPE_INTRON:
      gt_assert(!pw->processing_intron_with_1_base_left);
      gt_assert(!pw->processing_intron_with_2_bases_left);
      pw->gen_distance++;
      break;
    case EOP_TYPE_INTRON_WITH_1_BASE_LEFT:
      gt_assert(pw->proteineop);
      gt_assert(!pw->processing_intron_with_2_bases_left);
      if (!pw->processing_intron_with_1_base_left) {
        pw->processing_intron_with_1_base_left = true;
        pw->gen_distance++;
      }
      pw->gen_distance++;
      break;
    case EOP_TYPE_INTRON_WITH_2_BASES_LEFT:
      gt_assert(pw->proteineop);
      gt_assert(!pw->processing_intron_with_1_base_left);
      if (!pw->processing_intron_with_2_bases_left) {
        pw->processing_intron_with_2_bases_left = true;
        pw->gen_distance += 2;
        pw->ref_distance++;
      }
      pw->gen_distance++;
      break;
    case EOP_TYPE_MISMATCH_WITH_1_GAP:
      gt_assert(pw->proteineop);
      if (pw->processing_intron_with_1_base_left) {
        pw->processing_intron_with_1_base_left = false;
        pw->gen_distance++;
        pw->ref_distance++;
      }
      else if (pw->processing_intron_with_2_bases_left)
        pw->processing_intron_with_2_bases_left = false;
      else {
        pw->gen_distance += 2;
        pw->ref_distance++;
      }
      break;
    case EOP_TYPE_MISMATCH_WITH_2_GAPS:
      gt_assert(pw->proteineop);
      if (pw->processing_intron_with_1_base_left) {
        pw->processing_intron_with_1_base_left = false;
        pw->ref_distance++;
      }
      else {
        pw->gen_distance++;
        pw->ref_distance++;
      }
      break;
    case EOP_TYPE_DELETION_WITH_1_GAP:
      gt_assert(pw->proteineop);
      if (pw->processing_intron_with_1_base_left) {
        pw->processing_intron_with_1_base_left = false;
        pw->gen_distance++;
      }
      else
        pw->gen_distance += 2;
      break;
    case EOP_TYPE_DELETION_WITH_2_GAPS:
      gt_assert(pw->proteineop);
      if (!pw->processing_intron_with_1_base_left)
        pw->gen_distance++;
      break;
    default: gt_assert(0);
  }

  pw->last_eop_length--;

  if (pw->last_eop_length)
    pw->steps_in_current_eop++;
  else {
    pw->actual_eops++;
    pw->steps_in_current_eop = 0;
  }
}

void gth_path_walker_next(GthPathWalker *pw)
{
  gt_assert(pw && gth_path_walker_has_next(pw));
  if (!pw->last_eop_length) {
    pw->last_eop_type   = gt_editoperation_type(*pw->eopptr, pw->proteineop);
    pw->last_eop_length = gt_editoperation_length(*pw->eopptr, pw->proteineop);
    if (pw->forward)
      pw->eopptr--;
    else
      pw->eopptr++;
  }
  step(pw);
}

void gth_path_walker_try_steps(GthPathWalker *pw, unsigned long steps)
{
  unsigned long i;
  gt_assert(pw);
  for (i = 0; gth_path_walker_has_next(pw) && i < steps; i++)
    gth_path_walker_next(pw);
}

unsigned long gth_path_walker_eop_distance(const GthPathWalker *pw)
{
  gt_assert(pw);
  return pw->eop_distance;
}

unsigned long gth_path_walker_gen_distance(const GthPathWalker *pw)
{
  gt_assert(pw);
  return pw->gen_distance;
}

unsigned long gth_path_walker_ref_distance(const GthPathWalker *pw)
{
  gt_assert(pw);
  return pw->ref_distance;
}

void gth_path_walker_show(const GthPathWalker *pw, GtFile *outfp)
{
  gt_assert(pw);
  gt_file_xprintf(outfp, "GthPathWalker: orientation=%s, eop_dist=%lu, "
                     "gen_dist=%lu, ref_dist=%lu, actual=%lu, steps_cur=%d\n",
                     pw->forward ? "forward" : "reverse",
                     pw->eop_distance, pw->gen_distance, pw->ref_distance,
                     pw->actual_eops, pw->steps_in_current_eop);
}

unsigned long gth_path_walker_actual_eops(const GthPathWalker *pw)
{
  gt_assert(pw);
  return pw->actual_eops;
}

unsigned int gth_path_walker_steps_in_current_eop(const GthPathWalker *pw)
{
  gt_assert(pw);
  return pw->steps_in_current_eop;
}
