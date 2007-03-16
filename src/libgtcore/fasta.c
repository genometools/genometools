/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtcore/fasta.h>
#include <libgtcore/xansi.h>

void fasta_show_entry(const char *description, const char *sequence,
                      unsigned long sequence_length, unsigned long width)
{
  unsigned long i, current_length;
  xputchar(FASTA_SEPARATOR);
  if (description)
    xfputs(description, stdout);
  xputchar('\n');
  for (i = 0, current_length = 0; i < sequence_length; i++, current_length++) {
    if (width && current_length == width) {
      xputchar('\n');
      current_length = 0;
    }
    xputchar(sequence[i]);
  }
  xputchar('\n');
}
