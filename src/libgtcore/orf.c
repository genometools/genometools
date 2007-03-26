/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <libgtcore/orf.h>
#include <libgtcore/range.h>
#include <libgtcore/undef.h>

#define CODONLENGTH 3
#define START_AMINO 'M'
#define STOP_AMINO  '*'

void determine_ORFs(Array *ranges, unsigned int framenum,
                    const char *frame, unsigned long framelen, Env *env)
{
  unsigned long i;
  Range orf;
  assert(ranges && framenum < 3 && frame);
  orf.start = UNDEF_ULONG;
  for (i = 0; i < framelen; i++) {
    if (orf.start == UNDEF_ULONG) {
      if (frame[i] == START_AMINO)
        orf.start = i * CODONLENGTH + framenum;
    }
    else {
      if (frame[i] == STOP_AMINO) {
        orf.end = i * CODONLENGTH + framenum + 2;
        array_add(ranges, orf, env);
        orf.start = UNDEF_ULONG;
      }
    }
  }
}
