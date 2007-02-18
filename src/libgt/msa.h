/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef MSA_H
#define MSA_H

/* the multiple sequence alignment (MSA) class */
typedef struct MSA MSA;

MSA*          msa_new(const char *MSA_filename, Error*);
unsigned long msa_consensus_distance(const MSA*);
unsigned long msa_sum_of_pairwise_scores(const MSA*);
void          msa_show(MSA*); /* on stdout */
void          msa_free(MSA*);

#endif
