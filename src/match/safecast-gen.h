/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef SAFECAST_GEN_H
#define SAFECAST_GEN_H
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>

#ifdef S_SPLINT_S
#define PRINTLLUcast(X)  ((unsigned int) (X))
#else
#define PRINTLLUcast(X)  ((unsigned long long) (X))
#endif

#define DECLARESAFECASTFUNCTION(FROMTYPE,FROMTYPEALIAS,TOTYPE,TOTYPEALIAS)\
        static TOTYPE safecast_ ## FROMTYPEALIAS ## _ ## TOTYPEALIAS(\
                                                        const char *filename,\
                                                        int line,\
                                                        FROMTYPE from)\
        {\
          if (sizeof (FROMTYPE) > sizeof (TOTYPE))\
          {\
            if (from > (FROMTYPE) UINT32_MAX)\
            {\
              fprintf(stderr,"%s, %d: %llu cannot be stored in 32bit word",\
                        filename,line,PRINTLLUcast(from));\
              exit(EXIT_FAILURE); /* programming error */\
            }\
          }\
          return (TOTYPE) from;\
        }

#define CALLCASTFUNC(FROMTYPEALIAS,TOTYPEALIAS,VAL)\
        safecast_ ## FROMTYPEALIAS ## _ ## TOTYPEALIAS(__FILE__,__LINE__,VAL)

#endif
