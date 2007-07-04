/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2001 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtcore/tooldriver.h>
#include <tools/gt_skproto.h>

int main(int argc, char *argv[])
{
  return tooldriver(gt_skproto, argc, argv);
}
