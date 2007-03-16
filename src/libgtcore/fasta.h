/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FASTA_H
#define FASTA_H

#define FASTA_SEPARATOR '>'

/* show a fasta entry on stdout, if width is != 0 the sequence is formatted
   accordingly */
void fasta_show_entry(const char *description, const char *sequence,
                      unsigned long sequence_length, unsigned long width);

#endif
