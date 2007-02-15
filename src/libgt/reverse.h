/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef REVERSE_H
#define REVERSE_H

/* reverse 'dna_seq' of length 'seqlen' in place */
int reverse_complement(char *dna_seq, unsigned long seqlen, Error *err);

#endif
