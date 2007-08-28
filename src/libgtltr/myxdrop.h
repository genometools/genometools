/* David Ellinghaus */

#ifndef MYXDROP_H
#define MYXDROP_H

//#include "types.h"

typedef struct
{
  int mat,
      mis,
      ins,
      del,
      gcd;                    // greatest common divisor
} Arbitraryscores;             // \Typedef{ArbitraryScores}

typedef struct
{
  int mis,
      ins,
      del;
} Arbitrarydistances;             // \Typedef{ArbitraryScores}

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
                                     Xdropbest * xdropbest,\
	                             ArrayFrontvalue * fronts,\
	                             Uchar * useq,\
				     Uchar * vseq,\
	                             Sint ulen,\
				     Sint vlen,\
	                             Xdropscore xdropbelowscore)

#define EVALXDROPARBITSCORESLEFT\
      void evalxdroparbitscoresleft(Arbitraryscores *arbitscores,\
                                    Xdropbest * xdropbest,\
                                    ArrayFrontvalue * fronts,\
                                    Uchar * useq,\
				    Uchar * vseq,\
                                    Sint ulen,\
				    Sint vlen,\
                                    Xdropscore xdropbelowscore)
*/
#endif
