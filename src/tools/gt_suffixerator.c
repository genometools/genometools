/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

int gt_suffixerator(int argc, const char **argv, Env *env)
{
  env_error_check(env);
  return parseargsandcallsuffixerator(argc, argv, env);
}
