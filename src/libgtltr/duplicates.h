/*
  Copyright (c) 2007 David Ellinghaus <dellinghaus@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

void removeduplicates(ArrayLTRboundaries *arrayLTRboundaries);

void removeoverlapswithlowersimilarity(
  ArrayLTRboundaries *arrayLTRboundaries,
  bool nooverlapallowed);
