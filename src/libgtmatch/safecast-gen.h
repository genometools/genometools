/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

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
              exit(EXIT_FAILURE);\
            }\
          }\
          return (TOTYPE) from;\
        }

#define CALLCASTFUNC(FROMTYPEALIAS,TOTYPEALIAS,VAL)\
        safecast_ ## FROMTYPEALIAS ## _ ## TOTYPEALIAS(__FILE__,__LINE__,VAL)
