/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

#define SCAN_ALPHA_VALUE(NUM, CHAR_1, CHAR_2)                                  \
        if (sscanf(argv[NUM], "%d", &rval) != 1 || rval > 0) {                 \
          error("argument for alpha(%c,%c) must be non-positive integer",      \
                CHAR_1, CHAR_2);                                               \
        }                                                                      \
        scorematrix_set_score(energy_function, alpha_encode(dna_alpha, CHAR_1),\
                              alpha_encode(dna_alpha, CHAR_2), rval);          \
        scorematrix_set_score(energy_function, alpha_encode(dna_alpha, CHAR_2),\
                              alpha_encode(dna_alpha, CHAR_1), rval)

static int computeEentry(unsigned long i, unsigned long j, int **E,
                         char *rna_sequence, unsigned long rna_length,
                         ScoreMatrix *energy_function)
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
  alphavalue = scorematrix_get_score(energy_function, rna_sequence[i-1],
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
                           ScoreMatrix *energy_function)
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
        E[i][j] = computeEentry(i, j, E, rna_sequence, rna_length,
                                energy_function);
    }
  }
}

static void traceback(unsigned long i, unsigned long j, int **E,
                      char *rna_sequence, unsigned long rna_length,
                      ScoreMatrix *energy_function, FILE *fp)
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
      alphavalue = scorematrix_get_score(energy_function, rna_sequence[i-1],
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
                              ScoreMatrix *energy_function,
                              Alpha *dna_alpha, FILE *fp)
{
  unsigned long i;
  int **E;

  if (verbose) {
    fprintf(fp, "fold the following RNA sequence with Nussinov Algorithm:\n");
    for (i = 0; i < rna_length; i++) {
      xfputc(alpha_decode(dna_alpha, rna_sequence[i]), fp);
    }
    xfputc('\n', fp);
    fprintf(fp, "length of RNA sequence = %lu\n", rna_length);
    fprintf(fp, "using the following parameters:\n");
    fprintf(fp, "l_min =  %u\n", l_min);
    fprintf(fp, "alpha(G,C) = alpha(C,G) = %d\n",
            scorematrix_get_score(energy_function, alpha_encode(dna_alpha, 'G'),
                                  alpha_encode(dna_alpha, 'C')));
    fprintf(fp, "alpha(A,U) = alpha(U,A) = %d\n",
            scorematrix_get_score(energy_function, alpha_encode(dna_alpha, 'A'),
                                  alpha_encode(dna_alpha, 'U')));
    fprintf(fp, "alpha(G,U) = alpha(U,G) = %d\n",
            scorematrix_get_score(energy_function, alpha_encode(dna_alpha, 'U'),
                                  alpha_encode(dna_alpha, 'G')));
    fprintf(fp, "all other alpha values have been set to infinity\n");
  }

  /* alloc space for the matrix E */
  array2dim_calloc(E, rna_length+1, rna_length+1, int); /* XXX */

  /* compute matrix */
  compute_matrix(E, rna_sequence, rna_length, l_min, energy_function);

  /* perform backtracking to get optimal structure */
  if (verbose) {
    fprintf(fp, "score: %d\n", E[1][rna_length]);
    fprintf(fp, "result:\n");
  }
  traceback(1, rna_length, E, rna_sequence, rna_length, energy_function, fp);
  xfputc('\n', fp);

  /* free matrix E */
  array2dim_free(E);
}

int gt_nussinov_rna_fold(int argc, char *argv[], Error *err)
{
  unsigned long i, j, rna_length;
  unsigned int l_min;
  char *rna_sequence;
  int rval;
  Alpha *dna_alpha;
  ScoreMatrix *energy_function; /* alpha */
  error_check(err);

  /* check if correct arguments are given */
  if (argc !=6) {
    fprintf(stderr, "Usage: l_min alpha(G,C) alpha(A,U) alpha(G,U) "
                    "RNA_sequence\n");
    fprintf(stderr, "Fold the supplied RNA sequence with the Nussinov "
                    "algorithm.\n");
    exit(EXIT_FAILURE); /* XXX */
  }

  /* set DNA alphabet */
  dna_alpha = alpha_new_dna();
  energy_function = scorematrix_new(dna_alpha);

  /* init the energy function */
  for (i = 0; i < alpha_size(dna_alpha); i++) {
    for (j = 0; j < alpha_size(dna_alpha); j++) {
      scorematrix_set_score(energy_function, i, j, INT_MAX);
    }
  }

  /* save l_min value */
  if (sscanf(argv[1], "%d", &rval) != 1 || rval <= 0)
    error("argument for l_min must be positive integer");
  l_min = (unsigned int) rval;

  /* save alpha values */
  SCAN_ALPHA_VALUE(2, 'G', 'C');
  SCAN_ALPHA_VALUE(3, 'A', 'U');
  SCAN_ALPHA_VALUE(4, 'G', 'U');

  /* save RNA sequence */
  rna_length = strlen(argv[5]);
  rna_sequence = xstrdup(argv[5]);
  alpha_encode_seq(dna_alpha, rna_sequence, rna_sequence, rna_length);

  /* fold RNA sequence with Nussinov algorithm */
  nussinov_rna_fold(rna_sequence, rna_length, l_min, 1, energy_function,
                    dna_alpha, stdout);

  /* free */
  free(rna_sequence);
  scorematrix_free(energy_function);
  alpha_free(dna_alpha);

  return 0;
}
