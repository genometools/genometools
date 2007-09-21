/*
  Copyright (C) 2007 David Ellinghaus <dellinghaus@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef SEARCHFORLTRS_H
#define SEARCHFORLTRS_H

#include "libgtcore/env.h"
#include "libgtmatch/sarr-def.h"

#include "ltrharvest-opt.h"

int searchforLTRs(Sequentialsuffixarrayreader *ssar, LTRharvestoptions *lo,
                  const Seqpos *markpos, Env *env);

#endif
