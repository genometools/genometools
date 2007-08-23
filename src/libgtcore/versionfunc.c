/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdio.h>
#include "gt_build.h"
#include "gt_cc.h"
#include "gt_cflags.h"
#include "gt_version.h"

void versionfunc(const char *progname)
{
  printf("%s (GenomeTools) %s (%s)\n", progname, GT_VERSION, GT_BUILT);
  printf("Copyright (c) 2003-2007 Gordon Gremme, Stefan Kurtz, and "
         "CONTRIBUTORS\n");
  printf("Copyright (c) 2003-2007 Center for Bioinformatics, University of "
         "Hamburg\n");
  printf("See LICENSE file or http://genometools.org/license.html for license "
         "details.\n\n");
  printf("Used compiler: %s\n", GT_CC);
  printf("Compile flags: %s\n", GT_CFLAGS);
}

void showshortversion(const char *progname)
{
  printf("%s (GenomeTools) %s (%s)\n", progname, GT_VERSION, GT_BUILT);
}
