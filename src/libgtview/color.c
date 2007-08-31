/*
  Copyright (c) 2007 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "libgtview/color.h"

bool color_equals(Color c1, Color c2)
{
  return ((c1.red == c2.red) && (c1.green == c2.green) && (c1.blue == c2.blue));
}
