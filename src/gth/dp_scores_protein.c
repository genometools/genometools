/*
  Copyright (c) 2009-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#include "core/array2dim_api.h"
#include "core/array3dim.h"
#include "core/trans_table.h"
#include "core/unused_api.h"
#include "gth/align_common.h"
#include "gth/dp_scores_protein.h"

#define INDEL_PENALTY   -10.0
#define SCALEFACTOR     0.4

static GthFlt get_score(GtScoreMatrix *score_matrix,
                        GtAlphabet *score_matrix_alphabet,
                        unsigned char amino,
                        unsigned char origreferencechar)
{
  GthFlt rval = 0.0,
         scalefactor   = SCALEFACTOR,
         indel_penalty = INDEL_PENALTY;

  if (amino  == DASH || origreferencechar == DASH) {
    /* 1.) scaled INDEL_PENALTY for deletions from and insertions into genomic
       DNA of lengths 1, 2, or 3, irrespective of indel size */
    rval = scalefactor * indel_penalty;
  }
  else if (amino != WILDCARD && amino <= CHAR_MAX &&
           gt_alphabet_valid_input(score_matrix_alphabet, amino) &&
           origreferencechar <= CHAR_MAX &&
           gt_alphabet_valid_input(score_matrix_alphabet,
                                   origreferencechar)) {
    /* XXX: shorten this */
    if (amino == GT_STOP_AMINO) {
      /* 2.) (-)2*INDEL_PENALTY for matching/mismatching a stop codon */
      if (origreferencechar == GT_STOP_AMINO)
        rval = scalefactor * -2 * indel_penalty;
      else
        rval = scalefactor *  2 * indel_penalty;
    }
    else {
      /* 3.) amino acid substitution score */
      if (origreferencechar == GT_STOP_AMINO)
        rval = scalefactor *  2 * indel_penalty;
      else {
        GtUchar code1, code2;
        int wcidx;
        code1 = gt_alphabet_encode(score_matrix_alphabet, amino);
        code2 = gt_alphabet_encode(score_matrix_alphabet, origreferencechar);
        wcidx = gt_alphabet_size(score_matrix_alphabet) - 1;
        rval = scalefactor *
               gt_score_matrix_get_score(score_matrix,
                                         code1 == WILDCARD ? wcidx : code1,
                                         code2 == WILDCARD ? wcidx : code2);
      }
    }
  }
  /* 4.) else: neutral score in case of wild-card characters in the genomic DNA
   */

  return rval;
}

static GtUchar*** precompute_codon2amino(unsigned long translationtable)
{
  GtUchar ***codon2amino;
  GtTransTable *transtable;
  GtAlphabet *dna_alpha;
  int x, y, z, GT_UNUSED rval;
  char amino;
  gt_array3dim_malloc(codon2amino, 4, 4, 4);
  dna_alpha = gt_alphabet_new_dna();
  transtable = gt_trans_table_new(translationtable, NULL);
  /* XXX: the validity of the translation table has to be checked before */
  gt_assert(transtable);
  for (x = 0; x <= 3; x++) {
    for (y = 0; y <= 3; y++) {
      for (z = 0; z <= 3; z++) {
        char n1, n2, n3;
        n1 = gt_alphabet_decode(dna_alpha, x);
        n2 = gt_alphabet_decode(dna_alpha, y);
        n3 = gt_alphabet_decode(dna_alpha, z);
        rval = gt_trans_table_translate_codon(transtable, n1, n2, n3, &amino,
                                              NULL);
        gt_assert(!rval);
        codon2amino[x][y][z] = amino;
      }
    }
  }
  gt_trans_table_delete(transtable);
  gt_alphabet_delete(dna_alpha);
  return codon2amino;
}

static GthFlt** precompute_scores(GtScoreMatrix *score_matrix,
                                           GtAlphabet *score_matrix_alphabet)
{
  GthFlt **score;
  int x, y;
  gt_array2dim_malloc(score, UCHAR_MAX+1, UCHAR_MAX+1);
  for (x = 0; x < UCHAR_MAX+1; x++) {
    for (y = 0; y < UCHAR_MAX+1; y++)
      score[x][y] = get_score(score_matrix, score_matrix_alphabet, x, y);
  }
  return score;
}

GthDPScoresProtein* gth_dp_scores_protein_new(unsigned long translationtable,
                                              GtScoreMatrix *score_matrix,
                                              GtAlphabet *score_matrix_alphabet)
{
  GthDPScoresProtein *dp_scores_protein = gt_malloc(sizeof *dp_scores_protein);
  dp_scores_protein->codon2amino = precompute_codon2amino(translationtable);
  dp_scores_protein->score = precompute_scores(score_matrix,
                                               score_matrix_alphabet);
  return dp_scores_protein;
}

void gth_dp_scores_protein_delete(GthDPScoresProtein *dp_scores_protein)
{
  if (!dp_scores_protein) return;
  gt_array3dim_delete(dp_scores_protein->codon2amino);
  gt_array2dim_delete(dp_scores_protein->score);
  gt_free(dp_scores_protein);
}
