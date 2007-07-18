/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <limits.h>
#include "libgtcore/array2dim.h"
#include "libgtcore/bioseq.h"
#include <libgtext/msa.h>

#define GAPSYMBOL '-'

struct MSA {
  Bioseq *bs;
};

MSA* msa_new(const char *MSA_filename, Env *env)
{
  unsigned long i, firstseqlen;
  int had_err = 0;
  MSA *msa;
  env_error_check(env);
  msa = env_ma_malloc(env, sizeof (MSA));
  msa->bs = bioseq_new(MSA_filename, env);
  if (!msa->bs)
    had_err = -1;
  if (!had_err) {
    /* make sure that the MSA contains at least two sequences */
    if (bioseq_number_of_sequences(msa->bs) < 2) {
      env_error_set(env, "the MSA file '%s' contains less then 2 sequences",
                MSA_filename);
      had_err = -1;
    }
  }
  if (!had_err) {
    /* make sure that all sequences have the same length */
    firstseqlen = bioseq_get_sequence_length(msa->bs, 0);
    for (i = 1; i < bioseq_number_of_sequences(msa->bs); i++) {
      if (bioseq_get_sequence_length(msa->bs, i) != firstseqlen) {
        env_error_set(env, "length of sequence %lu in the MSA file '%s' "
                      "differs from the first", i, MSA_filename);
        had_err = -1;
        break;
      }
    }
  }
  if (had_err) {
    msa_delete(msa, env);
    return NULL;
  }
  return msa;
}

static char** get_msa_array(Bioseq *bs, Env *env)
{
  unsigned long i;
  char **msa;
  assert(bs);
  msa = env_ma_malloc(env,
                      sizeof (const char*) * bioseq_number_of_sequences(bs));
  for (i = 0; i < bioseq_number_of_sequences(bs); i++)
    msa[i] = (char*) bioseq_get_sequence(bs, i);
  return msa;
}

static unsigned long** get_count(char **msa, unsigned long number_of_seqs,
                                 unsigned long seqlen, Env *env)
{
  unsigned long col, seq, **count;
  assert(msa);
  array2dim_calloc(count, seqlen, UCHAR_MAX, unsigned long, env);
  for (seq = 0; seq < number_of_seqs; seq++) {
    for (col = 0; col < seqlen; col++)
      count[col][(int) msa[seq][col]]++;
  }
  return count;
}

static char* get_consensus(unsigned long **count, unsigned long seqlen,
                           Env *env)
{
  unsigned long col, c, max_count;
  char *consensus, consensus_char = GAPSYMBOL;
  assert(count);
  consensus = env_ma_malloc(env, sizeof (char) * seqlen);
  for (col = 0; col < seqlen; col++) {
    max_count = 0;
    for (c = 0; c < UCHAR_MAX; c++) {
      if (c != GAPSYMBOL && count[col][c] > max_count) {
        max_count = count[col][c];
        consensus_char = c;
      }
    }
    assert(consensus_char != GAPSYMBOL);
    consensus[col] = consensus_char;
  }
  return consensus;
}

unsigned long msa_consensus_distance(const MSA *msa, Env *env)
{
  unsigned long col, number_of_seqs, seqlen, **count, dist = 0;
  char **msa_array, *consensus;

  assert(msa);

  number_of_seqs = bioseq_number_of_sequences(msa->bs);
  seqlen = bioseq_get_sequence_length(msa->bs, 0);

  /* get the MSA in a convenient form */
  msa_array = get_msa_array(msa->bs, env);

  /* compute the character count array */
  count = get_count(msa_array, number_of_seqs, seqlen, env);

  /* compute the consensus from the count array */
  consensus = get_consensus(count, seqlen, env);

  /* compute the actual consensus distance */
  for (col = 0; col < seqlen; col++) {
    dist += number_of_seqs - count[col][(int) consensus[col]];
  }

  /* free */
  env_ma_free(consensus, env);
  array2dim_delete(count, env);
  env_ma_free(msa_array, env);

  return dist;
}

unsigned long msa_sum_of_pairwise_scores(const MSA *msa, Env *env)
{
  unsigned long i, j, col, number_of_seqs, seqlen, sum = 0;
  char **msa_array;

  assert(msa);

  number_of_seqs = bioseq_number_of_sequences(msa->bs);
  seqlen = bioseq_get_sequence_length(msa->bs, 0);

  /* get the MSA in a convenient form */
  msa_array = get_msa_array(msa->bs, env);

  /* compute the actual sum of pairwise scores */
  for (i = 0; i < number_of_seqs-1; i++) {
    for (j = i+1; j < number_of_seqs; j++) {
      for (col = 0; col < seqlen; col++)
        sum += (msa_array[i][col] == msa_array[j][col]) ? 0 : 1; /* delta */
    }
  }

  env_ma_free(msa_array, env);
  return sum;
}

void msa_show(MSA *msa)
{
  assert(msa);
  bioseq_show_as_fasta(msa->bs, 0);
}

void msa_delete(MSA *msa, Env *env)
{
  if (!msa) return;
  bioseq_delete(msa->bs, env);
  env_ma_free(msa, env);
}
