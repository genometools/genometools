/*
  Copyright (c) 2007 David Ellinghaus <d.ellinghaus@ikmb.uni-kiel.de>
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

#ifndef MYXDROP_H
#define MYXDROP_H

#include "libgtcore/arraydef.h"
#include "libgtcore/unused.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/seqpos-def.h"

typedef struct
{
  int mat,
      mis,
      ins,
      del,
      gcd;  /* greatest common divisor */
} Arbitraryscores;

typedef struct
{
  int mis,
      ins,
      del;
} Arbitrarydistances;

/*
  For each entry in the DP-matrix we store a single byte, and
  use the three rightmost bits to mark which edge in the edit distance
  graph to trace back.
*/

#define MYREPLACEMENTBIT   ((unsigned char) 1)          /* replacement */
#define MYDELETIONBIT      (((unsigned char) 1) << 1)   /* deletion */
#define MYINSERTIONBIT     (((unsigned char) 1) << 2)   /* insertion */
#define MYLENBIT           (((unsigned char) 1) << 3)
                                     /* length of identical substring */
#define MYMATCHBIT         (((unsigned char) 1) << 4)   /* match */
#define MYMISMATCHBIT      (((unsigned char) 1) << 5)   /* mismatch */

typedef struct
{
  int dptabrow;
  unsigned char dptabdirection; /* one of the bits REPLACEMENTBIT,
                                                   DELETIONBIT,
                                                   INSERTIONBIT */
} Myfrontvalue;

DECLAREARRAYSTRUCT(Myfrontvalue);

typedef struct
{
    unsigned int ivalue, jvalue;
    int score;
} Myxdropbest;

/* This is the type for the xdrop scores. */
typedef int Xdropscore;

DECLAREARRAYSTRUCT(Xdropscore);

int showmatrix ( ArrayMyfrontvalue * fronts,
  int distance,
  unsigned char * useq,
  unsigned char * vseq,
  int ulen,
  int vlen);

void calculatedistancesfromscores(Arbitraryscores *arbitscores,
    Arbitrarydistances *arbitdistances);

void calculateallowedMININFINITYINTgenerations(
   int *allowedMININFINITYINTgenerations,
   Arbitrarydistances *arbitdistances);

/*
   The following functions extend seeds to the right and to the left,
   respectively. xdropbest stores information about the best match
   found. useq is the first sequence position and vseq is the
   second sequence position. ulen and vlen are the
   remaining sequence length to align. If an alignment has score smaller than
   xdropbelowscore, then this alignment is not extended
   any more.
*/

void evalxdroparbitscoresright(Arbitraryscores *arbitscores,
                                     Myxdropbest * xdropbest,
                                     ArrayMyfrontvalue * fronts,
                                     const Encodedsequence *str_useq,
                                     const Encodedsequence *str_vseq,
                                     Seqpos useq,
                                     Seqpos vseq,
                                     int ulen,
                                     int vlen,
                                     Xdropscore xdropbelowscore,
                                     Error *err);

#define EVALXDROPARBITSCORESRIGHT\
      void evalxdroparbitscoresright(Arbitraryscores *arbitscores,\
                                     Myxdropbest * xdropbest,\
                                     ArrayMyfrontvalue * fronts,\
                                     const Encodedsequence *str_useq,\
                                     const Encodedsequence *str_vseq,\
                                     Seqpos useq,\
                                     Seqpos vseq,\
                                     int ulen,\
                                     int vlen,\
                                     Xdropscore xdropbelowscore,\
                                     UNUSED Error *err)

void evalxdroparbitscoresleft(Arbitraryscores * arbitscores,
                                    Myxdropbest * xdropbest,
                                    ArrayMyfrontvalue * fronts,
                                    const Encodedsequence *str_useq,
                                    const Encodedsequence *str_vseq,
                                    Seqpos useq,
                                    Seqpos vseq,
                                    int ulen,
                                    int vlen,
                                    Xdropscore xdropbelowscore,
                                    Error *err);

#define EVALXDROPARBITSCORESLEFT\
       void evalxdroparbitscoresleft(Arbitraryscores * arbitscores,\
                                    Myxdropbest * xdropbest,\
                                    ArrayMyfrontvalue * fronts,\
                                    const Encodedsequence *str_useq,\
                                    const Encodedsequence *str_vseq,\
                                    Seqpos useq,\
                                    Seqpos vseq,\
                                    int ulen,\
                                    int vlen,\
                                    Xdropscore xdropbelowscore,\
                                    UNUSED Error *err)
#endif
