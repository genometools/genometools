/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef XTMPFILE_H
#define XTMPFILE_H

#include <stdio.h>

#define XTMPFILE_TEMPLATE "/tmp/genometools.XXXXXXXXXX"

/* mkstemp(3) like function which returns a tmp file opened for writing */
FILE* xtmpfile(char *template);

#endif
