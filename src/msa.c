/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <limits.h>
#include "array2dim.h"
#include "bioseq.h"
#include "error.h"
#include "msa.h"
#include "undef.h"
#include "xansi.h"

#define GAPSYMBOL '-'

struct MSA {
  Bioseq *bs;
};

MSA* msa_new(const char *MSA_filename)
{
  unsigned long i, firstseqlen;
  MSA *msa;
  msa = xmalloc(sizeof(MSA));
  msa->bs = bioseq_new(MSA_filename);
  /* make sure that the MSA contains at least two sequences */
  if (bioseq_number_of_sequences(msa->bs) < 2)
    error("the MSA file '%s' contains less then 2 sequences", MSA_filename);
  /* make sure that all sequences have the same length */
  firstseqlen = bioseq_get_sequence_length(msa->bs, 0);
  for (i = 1; i < bioseq_number_of_sequences(msa->bs); i++) {
    if (bioseq_get_sequence_length(msa->bs, i) != firstseqlen) {
      error("length of sequence %lu in the MSA file '%s' differs from the "
            "first", i, MSA_filename);
    }
  }
  return msa;
}

static char** get_msa_array(Bioseq *bs)
{
  unsigned long i;
  char **msa;
  assert(bs);
  msa = xmalloc(sizeof(const char*) * bioseq_number_of_sequences(bs));
  for (i = 0; i < bioseq_number_of_sequences(bs); i++)
    msa[i] = (char*) bioseq_get_sequence(bs, i);
  return msa;
}

static unsigned long** get_count(char **msa, unsigned long number_of_seqs,
                                 unsigned long seqlen)
{
  unsigned long col, seq, **count;
  assert(msa);
  ARRAY2DIM_CALLOC(count, seqlen, UCHAR_MAX, unsigned long); 
  for (seq = 0; seq < number_of_seqs; seq++) { 
    for (col = 0; col < seqlen; col++)
      count[col][(int) msa[seq][col]]++;
  }
  return count;
}

static char* get_consensus(unsigned long **count, unsigned long seqlen)
{
  unsigned long col, c, max_count;
  char *consensus, consensus_char = GAPSYMBOL;
  assert(count);
  consensus = xmalloc(sizeof(char) * seqlen);
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

unsigned long msa_consensus_distance(const MSA *msa)
{
  unsigned long col, number_of_seqs, seqlen, **count, dist = 0;
  char **msa_array, *consensus;
  assert(msa);
  number_of_seqs = bioseq_number_of_sequences(msa->bs);
  seqlen = bioseq_get_sequence_length(msa->bs, 0);
  msa_array = get_msa_array(msa->bs);
  count = get_count(msa_array, number_of_seqs, seqlen);
  consensus = get_consensus(count, seqlen);
  for (col = 0; col < seqlen; col++) {
    dist += number_of_seqs - count[col][(int) consensus[col]];
  }
  free(consensus);
  ARRAY2DIM_FREE(count);
  free(msa_array);
  return dist;
}

unsigned long msa_sum_of_pairwise_scores(const MSA *msa)
{
  unsigned long i, j, col, number_of_seqs, seqlen, sum = 0;
  char **msa_array;
  assert(msa);
  number_of_seqs = bioseq_number_of_sequences(msa->bs);
  seqlen = bioseq_get_sequence_length(msa->bs, 0);
  msa_array = get_msa_array(msa->bs);
  for (i = 0; i < number_of_seqs-1; i++) {
    for (j = i+1; j < number_of_seqs; j++) {
      for (col = 0; col < seqlen; col++)
        sum += (msa_array[i][col] == msa_array[j][col]) ? 0 : 1; /* delta */
    }
  }
  free(msa_array);
  return sum;
}

void msa_show(MSA *msa)
{
  assert(msa);
  bioseq_show_as_fasta(msa->bs, 0);
}

void msa_free(MSA *msa)
{
  if (!msa);
  bioseq_free(msa->bs);
  free(msa);
}
