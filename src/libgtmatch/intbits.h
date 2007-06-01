#ifndef INTBITS_H
#define INTBITS_H
#include <limits.h>
#include "types.h"
#include "spacedef.h"

/*
  This file contains some definitions manipulating bitvectors represented
  by a \texttt{Uint}. In the comment lines we use $w$ for the word size
  and \texttt{\symbol{94}} for exponentiation of the previous character.
*/

#define INTWORDSIZE\
        (UintConst(1) << LOGWORDSIZE)     /* # of bits in Uint = w */
#define FIRSTBIT\
        (UintConst(1) << (INTWORDSIZE-1)) /* \(10^{w-1}\) */
#define ISBITSET(S,I)\
        (((S) << (I)) & FIRSTBIT)         /* is \(i\)th bit set? */
#define ITHBIT(I)\
        (FIRSTBIT >> (I))                 /* \(0^{i}10^{w-i-1}\)  */
#define SECONDBIT\
        (FIRSTBIT >> 1)                   /* \(010^{w-2}\) */
#define THIRDBIT\
        (FIRSTBIT >> 2)                   /* \(0010^{w-3}\) */
#define FOURTHBIT\
        (FIRSTBIT >> 3)                   /* \(00010^{w-4}\) */
#define FIFTHBIT\
        (FIRSTBIT >> 4)                   /* \(000010^{w-3}\) */
#define FIRSTTWOBITS\
        (UintConst(3) << (INTWORDSIZE-2)) /* \(11^{w-2}\) */
#define EXCEPTFIRSTBIT\
        (~FIRSTBIT)                       /* \(01^{w-1}\) */
#define EXCEPTFIRSTTWOBITS\
        (EXCEPTFIRSTBIT >> 1)             /* \(001^{w-2}\) */
#define EXCEPTFIRSTTHREEBITS\
        (EXCEPTFIRSTBIT >> 2)             /* \(0001^{w-3}\) */
#define EXCEPTFIRSTFOURBITS\
        (EXCEPTFIRSTBIT >> 3)             /* \(00001^{w-4}\) */
#define DIVWORDSIZE(I)\
        ((I) >> LOGWORDSIZE)              /* \((I) div w\) */
#define MODWORDSIZE(I)\
        ((I) & (INTWORDSIZE-1))           /* \((I) mod w\) */
#define MULWORDSIZE(I)\
        ((I) << LOGWORDSIZE)              /* \((I) * w\) */

/*
  The following defines the maximal value which can be represented by
  a given number of bits.
*/

#define MAXVALUEWITHBITS(BITNUM)       ((UintConst(1) << (BITNUM)) - 1)

/*
  The following defines the number of integers for a bitvector with N bits.
*/

#define NUMOFINTSFORBITS(N)\
        ((DIVWORDSIZE(N) == 0)\
           ? UintConst(1) \
           : (UintConst(1) + DIVWORDSIZE((N) - UintConst(1))))

/*
  The following macro allocates a bitarray of \texttt{N} bits. All bits
  are off.
*/

#define INITBITTABGENERIC(TAB,OLDTAB,NUMOFBITS)\
        {\
          Uint *tabptr, tabsize = NUMOFINTSFORBITS(NUMOFBITS);\
          ALLOCASSIGNSPACE(TAB,OLDTAB,Uint,tabsize);\
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
          Uint *tabptr, tabsize = NUMOFINTSFORBITS(N);\
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
