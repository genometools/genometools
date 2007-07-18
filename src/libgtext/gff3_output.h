/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GFF3_OUTPUT_H
#define GFF3_OUTPUT_H

#include "libgtext/genome_feature.h"

/* output the leading part of a genome feature in GFF3 format (i.e., the part
   up to the attributes) */
void gff3_output_leading(GenomeFeature *gf, GenFile*);

#endif
