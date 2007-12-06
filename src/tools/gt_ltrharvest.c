/*
  Copyright (c) 2007 David Ellinghaus <dellinghaus@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "libgtcore/error.h"
#include "libgtltr/ltrharvest-run.h"
#include "tools/gt_ltrharvest.h"

int gt_ltrharvest(int argc, const char **argv, Error *err)
{
  error_check(err);
  return parseargsandcallltrharvest(argc, argv, err);
}
