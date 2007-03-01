/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef SEQID2FILEOPTIONS_H
#define SEQID2FILEOPTIONS_H

#include "option.h"
#include "regionmapping.h"
#include "str.h"

/* add the options -seqfile and -regionmapping to the given option parser */
void           seqid2file_options(OptionParser*, Str *seqfile,
                                  Str *regionmapping, Env*);

RegionMapping* seqid2file_regionmapping_new(Str *seqfile, Str *regionmapping,
                                            Env *env);

#endif
