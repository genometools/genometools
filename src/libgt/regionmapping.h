/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef REGIONMAPPING_H
#define REGIONMAPPING_H

#include <libgt/str.h>

/* maps a sequence-region to a sequence file */
typedef struct RegionMapping RegionMapping;

RegionMapping* regionmapping_new_mapping(Str *mapping_filename, Env*);
RegionMapping* regionmapping_new_seqfile(Str *sequence_filename, Env*);
int            regionmapping_get_raw_sequence(RegionMapping*, const char**,
                                              Str *seqid, Env*);
int            regionmapping_get_raw_sequence_length(RegionMapping*,
                                                     unsigned long*,
                                                     Str *seqid, Env*);
void           regionmapping_delete(RegionMapping*, Env*);

#endif
