/*
  Copyright (C) 2007 David Ellinghaus <dellinghaus@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef OUTPUTGFF3_H
#define OUTPUTGFF3_H

int printgff3format(
  LTRharvestoptions *lo,
  Sequentialsuffixarrayreader *ssar,
  const Seqpos *markpos,
  Env *env);

#endif
