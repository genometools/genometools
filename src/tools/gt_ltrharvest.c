/*
  Copyright (c) 2007 David Ellinghaus <dellinghaus@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"
#include "libgtltr/ltrharvest-run.h"

int gt_ltrharvest(int argc, const char **argv, Env *env)
{
  env_error_check(env);
  return parseargsandcallltrharvest(argc, argv, env);
}
