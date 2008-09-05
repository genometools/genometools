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

#ifndef INTBITS_TAB_H
#define INTBITS_TAB_H
#include <string.h>
#include "intbits.h"
#include "spacedef.h"

#define DIVWORDSIZE(I)\
        ((I) >> LOGWORDSIZE)              /* \((I) div w\) */

#define MODWORDSIZE(I)\
        ((I) & (INTWORDSIZE-1))           /* \((I) mod w\) */

#define MULWORDSIZE(I)\
        ((I) << LOGWORDSIZE)              /* \((I) * w\) */

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
          size_t tabsize = NUMOFINTSFORBITS(NUMOFBITS);\
          ALLOCASSIGNSPACE(TAB,OLDTAB,Bitstring,tabsize);\
          (void) memset(TAB,0,sizeof(Bitstring) * tabsize);\
        }

#define INITBITTAB(TAB,N) INITBITTABGENERIC(TAB,NULL,N)

/*
  The following macro inititalizes a bitarray such that all bits
  are off.
*/

#define CLEARBITTAB(TAB,N)\
        (void) memset(TAB,0,sizeof(Bitstring) * NUMOFINTSFORBITS(N))

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

/*
  \texttt{BITNUM2WORD(TAB,I)} delivers the integer containing
  the \texttt{I}-th bit.
*/

#define BITNUM2WORD(TAB,I)  (TAB)[DIVWORDSIZE(I)]

#endif
