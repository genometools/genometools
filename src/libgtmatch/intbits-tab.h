/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef INTBITS_TAB_H
#define INTBITS_TAB_H
#include <limits.h>
#include "intbits.h"
#include "spacedef.h"

/*
  The following defines the number of integers for a bitvector with N bits.
*/

#define NUMOFINTSFORBITS(N)\
        ((DIVWORDSIZE(N) == 0)\
           ? (size_t) 1 \
           : ((size_t) 1 + (size_t) DIVWORDSIZE((N) - 1)))

/*
  The following macro allocates a bitarray of \texttt{N} bits. All bits
  are off.
*/

#define INITBITTABGENERIC(TAB,OLDTAB,NUMOFBITS)\
        {\
          Bitstring *tabptr;\
          size_t tabsize = NUMOFINTSFORBITS(NUMOFBITS);\
          ALLOCASSIGNSPACE(TAB,OLDTAB,Bitstring,tabsize);\
          for (tabptr = TAB; tabptr < (TAB) + tabsize; tabptr++)\
          {\
            *tabptr = 0;\
          }\
        }

#define INITBITTAB(TAB,N) INITBITTABGENERIC(TAB,NULL,N)

/*
  The following macro inititalizes a bitarray such that all bits
  are off.
*/

#define CLEARBITTAB(TAB,N)\
        {\
          Bitstring *tabptr;\
          size_t tabsize = NUMOFINTSFORBITS(N);\
          for (tabptr = TAB; tabptr < TAB + tabsize; tabptr++)\
          {\
            *tabptr = 0;\
          }\
        }

/*
  \texttt{SETIBIT(TAB,I)} sets the \texttt{I}-th bit in bitarray
  \texttt{TAB} to 1.
*/

#define SETIBIT(TAB,I)    (TAB)[DIVWORDSIZE(I)] |= ITHBIT(MODWORDSIZE(I))

/*
  \texttt{UNSSETIBIT(TAB,I)} sets the \texttt{I}-th bit in bitarray
  \texttt{TAB} to 0.
*/

#define UNSETIBIT(TAB,I)  (TAB)[DIVWORDSIZE(I)] &= ~(ITHBIT(MODWORDSIZE(I)))

/*
  \texttt{ISIBITSET(TAB,I)} checks if the \texttt{I}-th bit in bitarray
  \texttt{TAB} is 1.
*/

#define ISIBITSET(TAB,I)  ((TAB)[DIVWORDSIZE(I)] & ITHBIT(MODWORDSIZE(I)))

#endif
