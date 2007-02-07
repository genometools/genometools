/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef ENSURE_H
#define ENSURE_H

#include "error.h"

/* the ensure makro used for unit tests */
#define ensure(e)                                                      \
        (e) ?  0 : error("ensure \"%s\" failed: file \"%s\", line %d", \
                         #e, __FILE__, __LINE__);

#endif
