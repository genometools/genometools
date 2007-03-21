/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtcore/xtmpfile.h>
#include <libgtcore/xposix.h>

FILE* xtmpfile(char *template)
{
  FILE *tmpfp;
  int tmpfd;
  assert(template);
  tmpfd = xmkstemp(template);
  tmpfp = xfdopen(tmpfd, "w");
  return tmpfp;
}
