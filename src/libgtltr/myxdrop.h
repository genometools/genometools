/*
  Copyright (C) 2007 David Ellinghaus <dellinghaus@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef MYXDROP_H
#define MYXDROP_H

#include "libgtmatch/sarr-def.h"
#include "libgtmatch/seqpos-def.h"
#include "libgtmatch/arraydef.h"

typedef struct
{
  int mat,
      mis,
      ins,
      del,
      gcd;  // greatest common divisor
} Arbitraryscores;

typedef struct
{
  int mis,
      ins,
      del;
} Arbitrarydistances;

//#include "ltrharvest-opt.h"

/*
  For each entry in the DP-matrix we store a single byte, and
  use the three rightmost bits to mark which edge in the edit distance
  graph to trace back.
*/

#define MYREPLACEMENTBIT   ((unsigned char) 1)          // replacement
#define MYDELETIONBIT      (((unsigned char) 1) << 1)   // deletion
#define MYINSERTIONBIT     (((unsigned char) 1) << 2)   // insertion
#define MYLENBIT           (((unsigned char) 1) << 3)   // length of identical substring
#define MYMATCHBIT         (((unsigned char) 1) << 4)   // match
#define MYMISMATCHBIT      (((unsigned char) 1) << 5)   // mismatch

typedef struct
{
  int dptabrow;
  unsigned char dptabdirection; // one of the bits REPLACEMENTBIT,
                                //                 DELETIONBIT,
                                //                 INSERTIONBIT
} Myfrontvalue;

DECLAREARRAYSTRUCT(Myfrontvalue);

typedef struct
{
    unsigned int ivalue, jvalue;
    int score;
} Myxdropbest;

// This is the type for the xdrop scores.
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
   respectively. \texttt{bestmatch} stores information about the best match
   found. \texttt{useq} is the first sequence and \texttt{vseq} is the
   second sequence. \texttt{ulen} and \texttt{vlen} are the
   remaining sequence length to align. If an alignment has score smaller than
   \(M-\texttt{xdropbelowscore}\), then this alignment is not extended
   any more. \(M\) is the maximal score achieved so far.
*/
/*
#define EVALXDROPARBITSCORESRIGHT\
      void evalxdroparbitscoresright(Arbitraryscores *arbitscores,\
                                     Myxdropbest * xdropbest,\
	                             ArrayMyfrontvalue * fronts,\
	                             Seqpos * useq,\
				     Seqpos * vseq,\
	                             int ulen,\
				     int vlen,\
	                             Xdropscore xdropbelowscore)

#define EVALXDROPARBITSCORESLEFT\
      void evalxdroparbitscoresleft(Arbitraryscores *arbitscores,\
                                    Myxdropbest * xdropbest,\
                                    ArrayMyfrontvalue * fronts,\
                                    Seqpos * useq,\
				    Seqpos * vseq,\
                                    int ulen,\
				    int vlen,\
                                    Xdropscore xdropbelowscore)
*/

void evalxdroparbitscoresleft(Suffixarray *suffixarray,
                                    Arbitraryscores * arbitscores,
                                    Myxdropbest * xdropbest,
                                    ArrayMyfrontvalue * fronts,
                                    Seqpos useq,
				    Seqpos vseq,
                                    int ulen,
				    int vlen,
                                    Xdropscore xdropbelowscore,
				    Env *env);

void evalxdroparbitscoresright(Suffixarray *suffixarray,
                                     Arbitraryscores *arbitscores,
                                     Myxdropbest * xdropbest,
	                             ArrayMyfrontvalue * fronts,
	                             Seqpos useq,
				     Seqpos vseq,
	                             int ulen,
				     int vlen,
	                             Xdropscore xdropbelowscore,
				     Env *env);
#endif
