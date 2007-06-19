/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef TOOLDRIVER_H
#define TOOLDRIVER_H

#include <libgtcore/env.h>

/* The tool driver module allows to compile a tool into a separate binary. This
   is mostly useful for legacy applications like GenomeThreader.
   The tool driver creates an Env object, calls <tool>, and reports errors.
   See below for example code to create a separate binary for the eval tool.
   XXX: change example to reflect the real gth application
*/
int tooldriver(int(*tool)(int argc, const char **argv, Env*),
               int argc, char *argv[]);

#if 0

#include <libgtcore/tooldriver.h>
#include <tools/gt_gff3.h>

int main(int argc, char *argv[])
{
  return tooldriver(gt_gff3, argc, argv);
}

#endif

#endif
