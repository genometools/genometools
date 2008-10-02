/*
  Copyright (c) 2005-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2008 Center for Bioinformatics, University of Hamburg

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

#include <limits.h>
#include "core/array2dim.h"
#include "core/cstr.h"
#include "core/option.h"
#include "core/score_matrix.h"
#include "core/versionfunc.h"
#include "core/xansi.h"
#include "tools/gt_nussinov_rna_fold.h"

#define SCAN_ALPHA_VALUE(NUM, CHAR_1, CHAR_2)                                 \
        if (!had_err && (sscanf(argv[NUM], "%d", &rval) != 1 || rval > 0)) {  \
          gt_error_set(err, "argument for alpha(%c,%c) must be non-positive " \
                         "integer", CHAR_1, CHAR_2);                          \
          had_err = -1;                                                       \
        }                                                                     \
        if (!had_err) {                                                       \
          gt_score_matrix_set_score(energy_function,                          \
                                    gt_alpha_encode(dna_alpha, CHAR_1),       \
                                    gt_alpha_encode(dna_alpha, CHAR_2),       \
                                rval);                                        \
          gt_score_matrix_set_score(energy_function,                          \
                                    gt_alpha_encode(dna_alpha, CHAR_2),       \
                                    gt_alpha_encode(dna_alpha, CHAR_1),       \
                                    rval);                                    \
        }

static int computeEentry(unsigned long i, unsigned long j, int **E,
                         char *rna_sequence, GtScoreMatrix *energy_function)
{
  unsigned long k;
  int minvalue, value, alphavalue;

  /* 1. */
  minvalue = E[i][j-1];

  /* 2. */
  value =  E[i+1][j];
  if (value < minvalue)
    minvalue = value;

  /* 3. */
  alphavalue = gt_score_matrix_get_score(energy_function, rna_sequence[i-1],
                                      rna_sequence[j-1]);
  if (alphavalue != INT_MAX)
    value =  E[i+1][j-1] + alphavalue;
  else
    value = INT_MAX;
  if (value < minvalue)
    minvalue = value;

  /* 4. */
  for (k = i+1; k <= j-1; k++) {
    value =  E[i][k-1] + E[k][j];
    if (value < minvalue)
      minvalue = value;
  }

  return minvalue;
}

static void compute_matrix(int **E, char *rna_sequence,
                           unsigned long rna_length, unsigned int l_min,
                           GtScoreMatrix *energy_function)
{
  unsigned long i, j, l, n = rna_length;

  /* init matrix */
  for (i = 1; i <= rna_length; i++) {
    for (j = 1; j <= rna_length; j++) {
      if (i == j || i == j + 1)
        E[i][j] = 0; /* set all values in the antidiagonal to zero */
      else
        E[i][j] = INT_MAX;
    }
  }

  for (l = 1; l <= n - 1; l++) {
    for (i = 1; i <= n - l; i++) {
      j = i + l;
      if (j - i <= l_min)
        E[i][j] = 0;
      else
        E[i][j] = computeEentry(i, j, E, rna_sequence, energy_function);
    }
  }
}

static void traceback(unsigned long i, unsigned long j, int **E,
                      char *rna_sequence, unsigned long rna_length,
                      GtScoreMatrix *energy_function, FILE *fp)
{
  unsigned long k;
  int alphavalue;
  if (i < j) {
    if (E[i][j] == E[i+1][j]) {
      /* 1. */
      traceback(i+1, j, E, rna_sequence, rna_length, energy_function, fp);
    }
    else if (E[i][j] == E[i][j-1]) {
      /* 2. */
      traceback(i, j-1, E, rna_sequence, rna_length, energy_function, fp);
    }
    else {
      /* 3. */
      alphavalue = gt_score_matrix_get_score(energy_function, rna_sequence[i-1],
                                          rna_sequence[j-1]);
      if (alphavalue != INT_MAX && E[i][j] == E[i+1][j-1] + alphavalue) {
        fprintf(fp, "(%lu,%lu)", i, j);
        traceback(i+1, j-1, E, rna_sequence, rna_length, energy_function, fp);
      }
      else {
        for (k = i+1; k <= j-1; k++) {
          /* 4. */
          if (E[i][j] == E[i][k-1] + E[k][j]) {
            traceback(i, k-1, E, rna_sequence, rna_length, energy_function, fp);
            traceback(k, j, E, rna_sequence, rna_length, energy_function, fp);
          }
        }
      }
    }
  }
}

static void nussinov_rna_fold(char *rna_sequence, unsigned long rna_length,
                              unsigned int l_min, unsigned int verbose,
                              GtScoreMatrix *energy_function,
                              GtAlpha *dna_alpha, FILE *fp)
{
  unsigned long i;
  int **E;

  if (verbose) {
    fprintf(fp, "fold the following RNA sequence with Nussinov Algorithm:\n");
    for (i = 0; i < rna_length; i++) {
      gt_xfputc(gt_alpha_decode(dna_alpha, rna_sequence[i]), fp);
    }
    gt_xfputc('\n', fp);
    fprintf(fp, "length of RNA sequence = %lu\n", rna_length);
    fprintf(fp, "using the following parameters:\n");
    fprintf(fp, "l_min =  %u\n", l_min);
    fprintf(fp, "alpha(G,C) = alpha(C,G) = %d\n",
            gt_score_matrix_get_score(energy_function,
                                   gt_alpha_encode(dna_alpha, 'G'),
                                   gt_alpha_encode(dna_alpha, 'C')));
    fprintf(fp, "alpha(A,U) = alpha(U,A) = %d\n",
            gt_score_matrix_get_score(energy_function,
                                   gt_alpha_encode(dna_alpha, 'A'),
                                   gt_alpha_encode(dna_alpha, 'U')));
    fprintf(fp, "alpha(G,U) = alpha(U,G) = %d\n",
            gt_score_matrix_get_score(energy_function,
                                   gt_alpha_encode(dna_alpha, 'U'),
                                   gt_alpha_encode(dna_alpha, 'G')));
    fprintf(fp, "all other alpha values have been set to infinity\n");
  }

  /* alloc space for the matrix E */
  gt_array2dim_calloc(E, rna_length+1, rna_length+1); /* XXX */

  /* compute matrix */
  compute_matrix(E, rna_sequence, rna_length, l_min, energy_function);

  /* perform backtracking to get optimal structure */
  if (verbose) {
    fprintf(fp, "score: %d\n", E[1][rna_length]);
    fprintf(fp, "result:\n");
  }
  traceback(1, rna_length, E, rna_sequence, rna_length, energy_function, fp);
  gt_xfputc('\n', fp);

  /* free matrix E */
  gt_array2dim_delete(E);
}

static OPrval parse_options(int *parsed_args, int argc, const char **argv,
                            GtError *err)
{
  GtOptionParser *op;
  OPrval oprval;
  gt_error_check(err);
  op = gt_option_parser_new("l_min alpha(G,C) alpha(A,U) alpha(G,U) "
                            "RNA_sequence",
                            "Fold the supplied RNA sequence with the Nussinov "
                            "algorithm.");
  gt_option_parser_set_min_max_args(op, 5, 5);
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gt_versionfunc,
                                  err);
  gt_option_parser_delete(op);
  return oprval;
}

int gt_nussinov_rna_fold(int argc, const char **argv, GtError *err)
{
  unsigned long i, j, rna_length;
  unsigned int l_min = 0;
  char *rna_sequence = NULL;
  int parsed_args, rval, had_err = 0;
  GtAlpha *dna_alpha;
  GtScoreMatrix *energy_function; /* alpha */
  gt_error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  gt_assert(parsed_args == 1);

  /* set DNA alphabet */
  dna_alpha = gt_alpha_new_dna();
  energy_function = gt_score_matrix_new(dna_alpha);

  /* init the energy function */
  for (i = 0; i < gt_alpha_size(dna_alpha); i++) {
    for (j = 0; j < gt_alpha_size(dna_alpha); j++) {
      gt_score_matrix_set_score(energy_function, i, j, INT_MAX);
    }
  }

  /* save l_min value */
  if (sscanf(argv[1], "%d", &rval) != 1 || rval <= 0) {
    gt_error_set(err, "argument for l_min must be positive integer");
    had_err = -1;
  }
  else
    l_min = (unsigned int) rval;

  /* save alpha values */
  SCAN_ALPHA_VALUE(2, 'G', 'C');
  SCAN_ALPHA_VALUE(3, 'A', 'U');
  SCAN_ALPHA_VALUE(4, 'G', 'U');

  if (!had_err) {
    /* save RNA sequence */
    rna_length = strlen(argv[5]);
    rna_sequence = gt_cstr_dup(argv[5]);
    gt_alpha_encode_seq(dna_alpha, rna_sequence, rna_sequence, rna_length);

    /* fold RNA sequence with Nussinov algorithm */
    nussinov_rna_fold(rna_sequence, rna_length, l_min, 1, energy_function,
                      dna_alpha, stdout);
  }

  /* free */
  gt_free(rna_sequence);
  gt_score_matrix_delete(energy_function);
  gt_alpha_delete(dna_alpha);

  return had_err;
}
