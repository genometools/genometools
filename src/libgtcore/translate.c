/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <libgtcore/codon.h>
#include <libgtcore/translate.h>

void translate_dna(Str *protein, const char *dnaseq, unsigned long dnalen,
                   unsigned int frame, Env *env)
{
  const char *dnaptr;
  assert(protein && !str_length(protein) && dnaseq && frame < 3);
  /* translate the DNA in forward direction */
  for (dnaptr = dnaseq + frame; dnaptr < dnaseq + dnalen - 2; dnaptr += 3) {
    str_append_char(protein, codon2amino(*dnaptr, *(dnaptr+1), *(dnaptr+2)),
                    env);
  }
}
