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

#ifndef SEQPOS_DEF_H
#define SEQPOS_DEF_H

#include <stdbool.h>
#include <inttypes.h>

/*
  The following is the central case distinction to accomodate
  code for 32 bit integers and 64 bit integers.
*/

#ifdef S_SPLINT_S
#define PRINTSeqposcast(X)  ((unsigned long) (X))
#else
#define PRINTSeqposcast(X)  (X)
#endif

#ifdef _LP64

typedef uint64_t Seqpos;         /* \Typedef{Seqpos} */
#ifdef S_SPLINT_S
#define FormatSeqpos "%lu"
#else
#define FormatSeqpos   "%" PRIu64
#endif

#else

#ifdef S_SPLINT_S
#define FormatSeqpos "%lu"
#else
#define Formatuint64_t "%" PRIu64
#endif

#ifdef BIGSEQPOS
typedef uint64_t Seqpos;         /* \Typedef{Seqpos} */
#ifdef S_SPLINT_S
#define FormatSeqpos "%lu"
#else
#define FormatSeqpos "%" PRIu64
#endif
#else
typedef uint32_t Seqpos;         /* \Typedef{Seqpos} */
#ifdef S_SPLINT_S
#define FormatSeqpos "%lu"
#else
#define FormatSeqpos "%" PRIu32
#endif
#define Seqposequalsunsignedint
#endif /* BIGSEQPOS */

#endif /* _LP64 */

/*
  The following structure stores information about special characters.
*/

typedef struct
{
  Seqpos specialcharacters,      /* total number of special syms */
         specialranges,          /* number of ranges with special syms */
         realspecialranges,
         lengthofspecialprefix,  /* number of specials at start of sequence */
         lengthofspecialsuffix;  /* number of specials at end of sequence */
} Specialcharinfo;

typedef struct
{
  Seqpos position,
         value;
} Largelcpvalue;

typedef struct
{
  bool defined;
  Seqpos valueseqpos;
} DefinedSeqpos;

typedef struct
{
  Seqpos offset,
         left,
         right;
} Lcpinterval;

#endif
