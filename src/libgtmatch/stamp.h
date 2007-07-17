/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef STAMP_H
#define STAMP_H

#define STAMP\
        printf("STAMP(%d,%s)\n",__LINE__,__FILE__);\
        (void) fflush(stdout)

#ifndef STAMP
#define STAMP
#endif

#endif
