/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef ENSURE_H
#define ENSURE_H

#include "env.h"

/* the ensure macro used for unit tests */
#define ensure(has_err, e)                                                   \
        if (!has_err) {                                                      \
          if (!(e)) {                                                        \
            env_error_set(env, "ensure \"%s\" failed: file \"%s\", line %d", \
                          #e, __FILE__, __LINE__);                           \
            has_err = -1;                                                    \
          }                                                                  \
        }

#endif
