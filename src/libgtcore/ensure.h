/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef ENSURE_H
#define ENSURE_H

#include <libgtcore/env.h>

/* the ensure macro used for unit tests */
#define ensure(had_err, e)                                                   \
  do {                                                                       \
        if (!had_err) {                                                      \
          if (!(e)) {                                                        \
            env_error_set(env, "ensure \"%s\" failed: file \"%s\", line %d", \
                          #e, __FILE__, __LINE__);                           \
            had_err = -1;                                                    \
          }                                                                  \
        }                                                                    \
  } while (0)

#endif /* ENSURE_H */
