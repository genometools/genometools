/*
  Copyright (C) 2007 David Ellinghaus <dellinghaus@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <string.h>

/* Prof.Kurtz files*/
#include "libgtmatch/divmodmul.h"
#include "libgtmatch/chardef.h"

#include "repeattypes.h"
#include "myxdrop.h"
#include "minmax.h"

#define MINUSINFINITYINT ((int)integermin)
#define ACCESSTOFRONT(D,K) ((unsigned int) (D) * (D) + (D) + (K))

/*
 The following function shows the matrix of the calculated fronts.
 */
/* CAUTION: fronts, that run over the matrix boundaries are not shown in 
   the printed matrix. 
 */
int showmatrix ( ArrayMyfrontvalue * fronts,
  int distance,
  unsigned char * useq,
  unsigned char * vseq,
  int ulen,
  int vlen)
{
  int i, j, k, d = distance + 1, filled = 0;
  unsigned int l;
  int integermax, integermin;

  integermax = MAX (ulen, vlen), integermin = -integermax;

/*
  DEBUG1 (2,
          "Num of allocates Frontvalues = fronts->nextfreeFrontvalue: %lu\n",
          (Showuint) fronts->nextfreeFrontvalue);
*/
  printf( "ACHTUNG: Fronten, die ueber die boundaries gehen,"
                  " erscheinen nicht in Matrix.\n");
  printf("matrix:\n");
  printf("        ");
  printf("%-3c ", vseq[0]);
  // print vseq
  for (i = 1; i < vlen; i++)
  {
    printf("%-3c ", vseq[i]);
  }

  for (i = 0; i <= ulen; i++)
  {
    printf("\n");
    // print useq
    if (i != 0)
    {
      printf("%-3c ", useq[i - 1]);
    }
    else
    {
      printf("    ");
    }
    for (j = 0; j <= vlen; j++)
    {
      d = distance + 1;
      k = i - j;
      for (l = 0; l < fronts->nextfreeMyfrontvalue; l++)
      {
        if (MINUSINFINITYINT == fronts->spaceMyfrontvalue[l].dptabrow)
        {
          continue;
        }
        if (i == fronts->spaceMyfrontvalue[l].dptabrow)
        {
          for (d = 0; d <= distance; d++)
          {
            if ((k < -d) || (k > d))
            {
              continue;
            }
            if (l == ACCESSTOFRONT (d, i - j))
            {
              printf("%-3d ", d);
              l = fronts->nextfreeMyfrontvalue;
              filled++;
              break;
            }
          }
        }           //end if
      }         // end for l
      if (d > distance)
      {
        printf("U   ");
      }
    }              // end for j
  }            //end for i

  printf("\n%.2f percent of matrix filled\n",
           filled * 100.00 / ((ulen + 1) * (vlen + 1)));

  return 0;
}

/* If showmatrix funktion = "True", then
    (A)->nextfree##TYPE = POS+1; must be set:


 ** The other case, saves more space 

#define STOREINARRAYFRONTS(A,POS,TYPE,VAL)\
        CHECKARRAYSPACEMULTI(A,TYPE,POS+1)\
        (A)->space##TYPE[POS] = VAL;\
        (A)->nextfree##TYPE = POS+1;
*/
 
/* 
 POS+1 necessary, since at POS = 0 and nextfree##TYPE = 0
 nothing will be allocated in addition!
 */
#define STOREINARRAYFRONTS(A,POS,TYPE,VAL)\
        CHECKARRAYSPACEMULTI(A,TYPE,POS+1)\
        (A)->space##TYPE[POS] = VAL;

/*
 The following macro swaps two values.
 */
#define SWAPMAYBE\
        if(m < n)\
        {\
          r = m;\
          m = n;\
          n = r;\
        }
/*
 The following macro determines the greatest common divisor.
 */
#define GGT\
        do\
        {\
          r = m % n;\
          m = n;\
          n = r;\
        } while (r != 0)

/*
 The following function calculates the distance from the given scores.
 */
void calculatedistancesfromscores(Arbitraryscores *arbitscores, 
    Arbitrarydistances *arbitdistances)
{
  unsigned int m, n, r;
  int mat, mis, ins, del;
  
  // if mat is odd double all scores
  if(MOD2((unsigned int)arbitscores->mat))
  {    
    mat = arbitscores->mat  * (int)2; 
    mis = arbitscores->mis  * (int)2; 
    ins = arbitscores->ins  * (int)2; 
    del = arbitscores->del  * (int)2; 
    //DEBUG0(3, "scores has been doubled for calculation of distances.\n");
  }
  else
  {
    mat = arbitscores->mat;
    mis = arbitscores->mis;
    ins = arbitscores->ins;
    del = arbitscores->del;
  }

  m = (unsigned int)(mat - mis);
  //DEBUG1(2, "a = %lu\n", (Showuint) m);
  n = (unsigned int)(mat/2 - ins);
  //DEBUG1(2, "b = %lu\n", (Showuint) n);
  SWAPMAYBE  
  GGT;

  n = (unsigned int)(mat/2 - del);
  //DEBUG1(2, "c = %lu\n", (Showuint) n);
  SWAPMAYBE  
  GGT;
    
  arbitscores->gcd = (int) m;
  //DEBUG1(2, "g = gcd(a,b,c) = %ld\n", (Showsint) arbitscores->gcd);
  
  arbitdistances->mis  = (mat - mis) / arbitscores->gcd;
  arbitdistances->ins  = (mat/2 - ins)  / arbitscores->gcd;
  arbitdistances->del  = (mat/2 - del)  / arbitscores->gcd;
  
  /*
  DEBUG0(2, "distances:\n");
  DEBUG1(2, "mis = %ld\n", (Showsint) arbitdistances->mis);
  DEBUG1(2, "ins = %ld\n", (Showsint) arbitdistances->ins);
  DEBUG1(2, "del = %ld\n", (Showsint) arbitdistances->del);
  */
}

/*
 The following function calculates the maximal allowed number of
 generations with all front values equal minus infinity.
 */
void calculateallowedMININFINITYINTgenerations(
   int *allowedMININFINITYINTgenerations,
   Arbitrarydistances *arbitdistances)
{
  *allowedMININFINITYINTgenerations = MAX(arbitdistances->mis,
      			                  arbitdistances->ins);
  *allowedMININFINITYINTgenerations = MAX(*allowedMININFINITYINTgenerations,
      			                  arbitdistances->del);
  (*allowedMININFINITYINTgenerations)--;
}

/*
 The following macro checks for wildcard and separator symbols.
 */
#define COMPARESYMBOLSSEP(I,J)\
        USEQ(a,I);\
        if(a == (unsigned char) SEPARATOR)\
        {\
          ulen = I;\
          break;\
        }\
        VSEQ(b,J);\
        if(b == (unsigned char) SEPARATOR)\
        {\
          vlen = J;\
          break;\
        }\
        if(a != b || a == (unsigned char) WILDCARD)\
        {\
          break;\
        }

#ifdef XDROPDEF_H
  #undef MATCHSCORE
  #undef HALFMATCHSCORE
  #undef MISMATCHSCORE
  #undef S
  #undef SETDBACK
#endif

#define MATCHSCORE arbitscores->mat
#define HALFMATCHSCORE arbitscores->mat/2
#define GCD arbitscores->gcd
#define S(K,D)\
        ((K) * HALFMATCHSCORE - (D) * GCD)
#define SETDBACK(XDROPBELOWSCORE)\
        (XDROPBELOWSCORE + HALFMATCHSCORE) / GCD + 1

//#define EVALXDROPARBITSCORES EVALXDROPARBITSCORESRIGHT
//alt #define USEQ(A,I) A = useq[(I)]
//alt #define VSEQ(A,J) A = vseq[(J)]
#define USEQ(A,I) A = getencodedchar(suffixarray->encseq,\
                                     useq+(Seqpos)(I),\
                                     Forwardmode)
#define VSEQ(A,J) A = getencodedchar(suffixarray->encseq,\
                                     vseq+(Seqpos)(J),\
				     Forwardmode)

//#include "myxdrop.gen"
void evalxdroparbitscoresright(Suffixarray *suffixarray,
                                     Arbitraryscores *arbitscores,
                                     Myxdropbest * xdropbest,
	                             ArrayMyfrontvalue * fronts,
	                             Seqpos useq,
				     Seqpos vseq,
	                             int ulen,
				     int vlen,
	                             Xdropscore xdropbelowscore,
				     Env *env)
{
  int integermax, integermin, i , // index i 
    j = 0,                            // index j
    k = 0,                            // diagonal
    end_k = ulen - vlen,              // diagonal of endpoint (ulen, vlen)
    lbound = 0,                       // diagonal lower bound
    ubound = 0,                       // diagonal upper bound
    lboundtmp = 0, 
    uboundtmp = 0, 
    tmprow;                           // temporary row index
  unsigned char a, b;
  Myfrontvalue tmpfront;
  Arbitrarydistances arbitdistances;
  int d = 0,                         // distance
      d_pre = 0,                     // distance previous
      dback;
  Xdropscore bigt_tmp;                // best score T' seen already
  ArrayXdropscore big_t;              // array T[d] of best score T for each d
  bool alwaysMININFINITYINT = true;
  int allowedMININFINITYINTgenerations = 0,
      currentMININFINITYINTgeneration = 0; 
  
  //DEBUG0(1, "starting xdrop-alignment-extension...\n");
  
  calculatedistancesfromscores(arbitscores, &arbitdistances);
  calculateallowedMININFINITYINTgenerations(
     &allowedMININFINITYINTgenerations,
     &arbitdistances);
  
  //DEBUG1(2, "allowedMININFINITYINTgenerations = %ld\n", 
  //    (Showsint) allowedMININFINITYINTgenerations);

  dback = SETDBACK(xdropbelowscore);
  //DEBUG1(2, "dback = %ld\n", (Showsint) dback);
  
  INITARRAY (&big_t, Xdropscore);
  integermax = MAX (ulen, vlen);
  integermin = -integermax;

  // "phase" 0
  for(i = 0; i < MIN(ulen,vlen); i++)
  {
    COMPARESYMBOLSSEP(i,i);
  }

  /* alignment already finished */
  if(i >= ulen || i >= vlen )
  {
    lbound =  1;
    ubound = -1;
  }

  /*
  DEBUG2 (3, "\nulen = %ld\nvlen = %ld\n",
          (Showsint) ulen, (Showsint) vlen);
  DEBUG5 (3, "nach phase %ld:   d = %ld   k = %ld   i = %ld   j = %ld\n",
          (Showsint) d, (Showsint) d, (Showsint) k, (Showsint) i,
          (Showsint) i);
  */

  tmpfront.dptabrow = i;
  tmpfront.dptabdirection = (unsigned char) 0;   // no predecessor

  STOREINARRAYFRONTS (fronts, 0, Myfrontvalue, tmpfront);
  xdropbest->score = bigt_tmp = S (i + i, 0);
  xdropbest->ivalue = xdropbest->jvalue = (unsigned int) i;

  /*
  DEBUG4 (3,
          "nach phase %ld:   best score so far = T' = T[0] ="
	  "S'(%ld,%ld) = %ld\n\n",
          (Showsint) d, (Showsint) i + i, (Showsint) d,
          (Showsint) bigt_tmp);
  */

  STOREINARRAY (&big_t, Xdropscore, 10, bigt_tmp);

  // "phase" d > 0
  while (lbound <= ubound)
  {
    d++;
    d_pre = d - dback; 

    // calculate fronts
    for (k = lbound - 1; k <= ubound + 1; k++)
    {
        
      i = MINUSINFINITYINT; 
      // case 1 : DELETION-EDGE
      if ( (lbound < k) && 
	   (d - arbitdistances.del >= 0) && 
	   (( -(d - arbitdistances.del) <= k - 1) && (k - 1 <= d - arbitdistances.del)) 
	 )
      {
        i =
          fronts->spaceMyfrontvalue[ACCESSTOFRONT(d - arbitdistances.del, k - 1)].dptabrow +
          1;
        tmpfront.dptabdirection = MYDELETIONBIT;
      }
      // case 2: REPLACEMENT-EDGE
      if ( (lbound <= k) && (k <= ubound) && (d - arbitdistances.mis >= 0) &&
	   (( -(d - arbitdistances.mis) <= k) && (k <= d - arbitdistances.mis)) 
	 )
      {
        // test, if case 1 has happened.
        if (tmpfront.dptabdirection & MYDELETIONBIT)
        {
          tmprow =
            fronts->spaceMyfrontvalue[ACCESSTOFRONT(d - arbitdistances.mis, k)].dptabrow + 1;
          if (tmprow > i)
          {
            i = tmprow;
            tmpfront.dptabdirection = MYREPLACEMENTBIT;
          }
        }
        else
        {
          i =
            fronts->spaceMyfrontvalue[ACCESSTOFRONT(d - arbitdistances.mis, k)].dptabrow + 1;
          tmpfront.dptabdirection = MYREPLACEMENTBIT;
        }

      }
      // case 3: INSERTION-EDGE
      if ( (k < ubound) && (d - arbitdistances.ins >= 0) &&
	   (( -(d - arbitdistances.ins) <= k + 1) && (k + 1 <= d - arbitdistances.ins)) 
	 )
      {
        if ((tmpfront.dptabdirection & MYDELETIONBIT)
            || (tmpfront.dptabdirection & MYREPLACEMENTBIT))
        {
          tmprow =
            fronts->spaceMyfrontvalue[ACCESSTOFRONT (d - arbitdistances.ins, k + 1)].dptabrow;
          if (tmprow > i)
          {
            i = tmprow;
            tmpfront.dptabdirection = MYINSERTIONBIT;
          }
        }
        else
        {
          i =
            fronts->spaceMyfrontvalue[ACCESSTOFRONT (d - arbitdistances.ins, k + 1)].dptabrow;
          tmpfront.dptabdirection = MYINSERTIONBIT;
        }
      }


      j = i - k;

      // if i = MINUSINFINITYINY or MINUSINFINITYINY + 1
      if( i < 0)
      {
	if(tmpfront.dptabdirection == (unsigned char)0)
	{
          alwaysMININFINITYINT = false;
	}
        tmpfront.dptabrow = MINUSINFINITYINT;
        /*
	  DEBUG5 (3,
                "nach phase %ld:   d = %ld   k = %ld   (i,j) = (%ld,%ld)\n",
                (Showsint) d, (Showsint) d, (Showsint) k, 
		(Showsint) i, (Showsint) i-k);
        */
      }
      // alignment score smaller than T - X
      else if( d_pre > 0 && 
      	      (S (i + j, d) < big_t.spaceXdropscore[d_pre] - xdropbelowscore)
	     )
      {
        tmpfront.dptabrow = MINUSINFINITYINT;
        /*
	  DEBUG5 (3,
                "nach phase %ld:   d = %ld   k = %ld   (i,j) = (%ld,%ld)\n",
                (Showsint) d, (Showsint) d, (Showsint) k, 
		(Showsint) i, (Showsint) i-k);
        */
      }
      else if( ( k <= -d || k >= d ) /* not correct boundaries for 
	                                ACCESTOFRONT(d-1,k) */
               ||
	       ( (fronts->spaceMyfrontvalue[ACCESSTOFRONT(d - 1, k)].dptabrow < i)
	         &&
	         (i <= MIN(ulen,vlen + k))
	       ) 
	     )
      {
        while(i < ulen && j < vlen) 
        {
	  COMPARESYMBOLSSEP(i,j);
          i++;
          j++;
        }
      
        alwaysMININFINITYINT = false;
	tmpfront.dptabrow = i;
     
	// MAX(bigt_tmp, S(i+j,d));
        if (S (i + j, d) > bigt_tmp)
        {
          xdropbest->score = bigt_tmp = S (i + j, d);
          xdropbest->ivalue = (unsigned int) i;
          xdropbest->jvalue = (unsigned int) j;
        }

        /*
	DEBUG5 (3,
                "nach phase %ld:   d = %ld   k = %ld  (i,j) = (%ld,%ld)\n",
                (Showsint) d, (Showsint) d, (Showsint) k, (Showsint) i,
                (Showsint) j);
        */
      }
       
      else 
      {
        alwaysMININFINITYINT = false;
        tmpfront.dptabrow = 
	  fronts->spaceMyfrontvalue[ACCESSTOFRONT(d - 1, k)].dptabrow;
        /*
        DEBUG5 (3,
                "nach phase %ld:   d = %ld   k = %ld   (i,j) = (%ld,%ld)\n",
                (Showsint) d, (Showsint) d ,(Showsint) k, 
		(Showsint) tmpfront.dptabrow, (Showsint) tmpfront.dptabrow-k);
        */
      }

 
      STOREINARRAYFRONTS (fronts, ACCESSTOFRONT (d, k), Myfrontvalue,
                            tmpfront);

      /*
      DEBUG3(3, "fronts(%ld,%ld) = rowindex = %ld\n", 
	  (Showsint) d,
	  (Showsint) k,
	  (Showsint) fronts->spaceMyfrontvalue[ACCESSTOFRONT(d,k)].dptabrow);
      */
      
      // delete value for test for 
      // INSERTIONBIT/REPLACEMENTBIT/DELETIONBIT above
      tmpfront.dptabdirection = (unsigned char) 0;


    }   // end for loop

    // if all front values are MINUSINFINITYINT, 
    //   aligment prematurely finished if allowedMININFINITYINTgenerations 
    //   exceeded
    //  (full front has already ended at d - currentMININFINITYINTgeneration).
    if (alwaysMININFINITYINT)
    {
      currentMININFINITYINTgeneration++;
      if(currentMININFINITYINTgeneration > allowedMININFINITYINTgenerations)
      {
	d = d - currentMININFINITYINTgeneration;
        /*
	DEBUG1(1,
	    "Alignment dropped prematurely because of threshold "
	    "at distance = "
	    "%ld\n",
	    (Showsint) d);
        */
	break;
      }
    }
    else
    {
      currentMININFINITYINTgeneration = 0;
      alwaysMININFINITYINT = true;
    }
    
    /*
    DEBUG4 (3,
            "nach phase %ld:   best score so far = T' = %ld, (i,j) = "
	    "(%lu,%lu)\n\n",
            (Showsint) d, (Showsint) bigt_tmp,
            (Showuint) xdropbest->ivalue, (Showuint) xdropbest->jvalue);
    */
    STOREINARRAY (&big_t, Xdropscore, 10, bigt_tmp);

    // fill out of bounds values of MINUSINFINITYINT
    // needed for showmatrix function
    for (k = -d; k < lbound - 1; k++)
    {
      tmpfront.dptabrow = MINUSINFINITYINT;
      STOREINARRAYFRONTS (fronts, ACCESSTOFRONT (d, k), Myfrontvalue,
                          tmpfront);
    }
    for (k = ubound + 2; k <= d; k++)
    {
      tmpfront.dptabrow = MINUSINFINITYINT;
      STOREINARRAYFRONTS (fronts, ACCESSTOFRONT (d, k), Myfrontvalue,
                          tmpfront);
    }

    // alignment finished 
    if ((-d <= end_k && end_k <= d) &&
        fronts->spaceMyfrontvalue[ACCESSTOFRONT (d, end_k)].dptabrow == ulen)
    {
      /*
      DEBUG3 (2,
              "nach phase %ld: endpoint (%ld,%ld) reached. "
	      "alignment finished.\n",
              (Showsint) d, (Showsint) ulen, (Showsint) vlen);
      */
      break;
    }

    // pruning lower bound
    // lbound may decrease by one or increase/stays the same
    for (k = lbound - 1; k <= ubound + 1; k++)
    {
      if (fronts->spaceMyfrontvalue[ACCESSTOFRONT (d, k)].dptabrow >
          MINUSINFINITYINT)
      {
        lbound = k;
        break;
      }
    }

    //pruning upper bound
    // ubound may increase by one or decrease/stays the same
    for (k = ubound + 1; k >= lbound - 1; k--)
    {
      if (fronts->spaceMyfrontvalue[ACCESSTOFRONT (d, k)].dptabrow >
          MINUSINFINITYINT)
      {
        ubound = k;
        break;
      }
    }

    //handling boundaries lower bound
    lboundtmp = lbound;
    for (k = 0; k >= lbound; k--)
    {
      if (fronts->spaceMyfrontvalue[ACCESSTOFRONT (d, k)].dptabrow == vlen + k)
      {
        // debugging 12.11.06, vorher:
        //lboundtmp = k + 2; //bei Zhang so
        lboundtmp = k;
        break;
      }
    }

    //handling boundaries upper bound
    uboundtmp = ubound;
    for (k = 0; k <= ubound; k++)
    {
      if (fronts->spaceMyfrontvalue[ACCESSTOFRONT (d, k)].dptabrow == ulen)
      {
        // debugging 12.11.06, vorher:
        //uboundtmp = k - 2; //bei Zhang so
        uboundtmp = k;
        break;
      }
    }

    lbound = MAX (lbound, lboundtmp);
    ubound = MIN (ubound, uboundtmp);
  }

  /*
  DEBUG0(1, "Alignment results:\n");
  DEBUG1(1, "distance = %ld\n", (Showsint) d);
  DEBUG2(1,
           "best (prefix-) alignment score from (i,j) = (0,0) to "
	   "(%lu,%lu):\n",
           (Showuint) xdropbest->ivalue, (Showuint) xdropbest->jvalue);
  DEBUG1(1, "score T' = %ld\n", (Showsint) bigt_tmp);
  */

  FREEARRAY (&big_t, Xdropscore);
}

#undef EVALXDROPARBITSCORES
#undef USEQ
#undef VSEQ

/*
  Now we redefine the macros to compute the left to right
  we use macros to abstract from the differences. The 
  following 4 macros are defined for the right to left extension beginning
  with the last character of the strings under consideration.
*/

/* CAUTION: special case for xdrop in LTRharvest,
   see following USEQ, VSEQ Macros */

//#define EVALXDROPARBITSCORES EVALXDROPARBITSCORESLEFT
//#define USEQ(A,I) A = *(useq-1-(I))
//#define VSEQ(A,J) A = *(vseq-1-(J))
#define USEQ(A,I) A = getencodedchar(suffixarray->encseq,\
                                     useq-1-(I),\
				     Forwardmode)
#define VSEQ(A,J) A = getencodedchar(suffixarray->encseq,\
                                     vseq-1-(J),\
				     Forwardmode)

void evalxdroparbitscoresleft(Suffixarray *suffixarray,
                                     Arbitraryscores *arbitscores,
                                     Myxdropbest * xdropbest,
	                             ArrayMyfrontvalue * fronts,
	                             Seqpos useq,
				     Seqpos vseq,
	                             int ulen,
				     int vlen,
	                             Xdropscore xdropbelowscore,
				     Env *env)
{
  int integermax, integermin, i , // index i 
    j = 0,                            // index j
    k = 0,                            // diagonal
    end_k = ulen - vlen,              // diagonal of endpoint (ulen, vlen)
    lbound = 0,                       // diagonal lower bound
    ubound = 0,                       // diagonal upper bound
    lboundtmp = 0, 
    uboundtmp = 0, 
    tmprow;                           // temporary row index
  unsigned char a, b;
  Myfrontvalue tmpfront;
  Arbitrarydistances arbitdistances;
  int d = 0,                         // distance
      d_pre = 0,                     // distance previous
      dback;
  Xdropscore bigt_tmp;                // best score T' seen already
  ArrayXdropscore big_t;              // array T[d] of best score T for each d
  bool alwaysMININFINITYINT = true;
  int allowedMININFINITYINTgenerations = 0,
      currentMININFINITYINTgeneration = 0; 
  
  //DEBUG0(1, "starting xdrop-alignment-extension...\n");
  
  calculatedistancesfromscores(arbitscores, &arbitdistances);
  calculateallowedMININFINITYINTgenerations(
     &allowedMININFINITYINTgenerations,
     &arbitdistances);
  
  //DEBUG1(2, "allowedMININFINITYINTgenerations = %ld\n", 
  //    (Showsint) allowedMININFINITYINTgenerations);

  dback = SETDBACK(xdropbelowscore);
  //DEBUG1(2, "dback = %ld\n", (Showsint) dback);
  
  INITARRAY (&big_t, Xdropscore);
  integermax = MAX (ulen, vlen);
  integermin = -integermax;

  // "phase" 0
  for(i = 0; i < MIN(ulen,vlen); i++)
  {
    COMPARESYMBOLSSEP(i,i);
  }

  /* alignment already finished */
  if(i >= ulen || i >= vlen )
  {
    lbound =  1;
    ubound = -1;
  }

  /*
  DEBUG2 (3, "\nulen = %ld\nvlen = %ld\n",
          (Showsint) ulen, (Showsint) vlen);
  DEBUG5 (3, "nach phase %ld:   d = %ld   k = %ld   i = %ld   j = %ld\n",
          (Showsint) d, (Showsint) d, (Showsint) k, (Showsint) i,
          (Showsint) i);
  */

  tmpfront.dptabrow = i;
  tmpfront.dptabdirection = (unsigned char) 0;   // no predecessor

  STOREINARRAYFRONTS (fronts, 0, Myfrontvalue, tmpfront);
  xdropbest->score = bigt_tmp = S (i + i, 0);
  xdropbest->ivalue = xdropbest->jvalue = (unsigned int) i;

  /*
  DEBUG4 (3,
          "nach phase %ld:   best score so far = T' = T[0] ="
	  "S'(%ld,%ld) = %ld\n\n",
          (Showsint) d, (Showsint) i + i, (Showsint) d,
          (Showsint) bigt_tmp);
  */

  STOREINARRAY (&big_t, Xdropscore, 10, bigt_tmp);

  // "phase" d > 0
  while (lbound <= ubound)
  {
    d++;
    d_pre = d - dback; 

    // calculate fronts
    for (k = lbound - 1; k <= ubound + 1; k++)
    {
        
      i = MINUSINFINITYINT; 
      // case 1 : DELETION-EDGE
      if ( (lbound < k) && 
	   (d - arbitdistances.del >= 0) && 
	   (( -(d - arbitdistances.del) <= k - 1) && (k - 1 <= d - arbitdistances.del)) 
	 )
      {
        i =
          fronts->spaceMyfrontvalue[ACCESSTOFRONT(d - arbitdistances.del, k - 1)].dptabrow +
          1;
        tmpfront.dptabdirection = MYDELETIONBIT;
      }
      // case 2: REPLACEMENT-EDGE
      if ( (lbound <= k) && (k <= ubound) && (d - arbitdistances.mis >= 0) &&
	   (( -(d - arbitdistances.mis) <= k) && (k <= d - arbitdistances.mis)) 
	 )
      {
        // test, if case 1 has happened.
        if (tmpfront.dptabdirection & MYDELETIONBIT)
        {
          tmprow =
            fronts->spaceMyfrontvalue[ACCESSTOFRONT(d - arbitdistances.mis, k)].dptabrow + 1;
          if (tmprow > i)
          {
            i = tmprow;
            tmpfront.dptabdirection = MYREPLACEMENTBIT;
          }
        }
        else
        {
          i =
            fronts->spaceMyfrontvalue[ACCESSTOFRONT(d - arbitdistances.mis, k)].dptabrow + 1;
          tmpfront.dptabdirection = MYREPLACEMENTBIT;
        }

      }
      // case 3: INSERTION-EDGE
      if ( (k < ubound) && (d - arbitdistances.ins >= 0) &&
	   (( -(d - arbitdistances.ins) <= k + 1) && (k + 1 <= d - arbitdistances.ins)) 
	 )
      {
        if ((tmpfront.dptabdirection & MYDELETIONBIT)
            || (tmpfront.dptabdirection & MYREPLACEMENTBIT))
        {
          tmprow =
            fronts->spaceMyfrontvalue[ACCESSTOFRONT (d - arbitdistances.ins, k + 1)].dptabrow;
          if (tmprow > i)
          {
            i = tmprow;
            tmpfront.dptabdirection = MYINSERTIONBIT;
          }
        }
        else
        {
          i =
            fronts->spaceMyfrontvalue[ACCESSTOFRONT (d - arbitdistances.ins, k + 1)].dptabrow;
          tmpfront.dptabdirection = MYINSERTIONBIT;
        }
      }


      j = i - k;

      // if i = MINUSINFINITYINY or MINUSINFINITYINY + 1
      if( i < 0)
      {
	if(tmpfront.dptabdirection == (unsigned char)0)
	{
          alwaysMININFINITYINT = false;
	}
        tmpfront.dptabrow = MINUSINFINITYINT;
        /*
	  DEBUG5 (3,
                "nach phase %ld:   d = %ld   k = %ld   (i,j) = (%ld,%ld)\n",
                (Showsint) d, (Showsint) d, (Showsint) k, 
		(Showsint) i, (Showsint) i-k);
        */
      }
      // alignment score smaller than T - X
      else if( d_pre > 0 && 
      	      (S (i + j, d) < big_t.spaceXdropscore[d_pre] - xdropbelowscore)
	     )
      {
        tmpfront.dptabrow = MINUSINFINITYINT;
        /*
	  DEBUG5 (3,
                "nach phase %ld:   d = %ld   k = %ld   (i,j) = (%ld,%ld)\n",
                (Showsint) d, (Showsint) d, (Showsint) k, 
		(Showsint) i, (Showsint) i-k);
        */
      }
      else if( ( k <= -d || k >= d ) /* not correct boundaries for 
	                                ACCESTOFRONT(d-1,k) */
               ||
	       ( (fronts->spaceMyfrontvalue[ACCESSTOFRONT(d - 1, k)].dptabrow < i)
	         &&
	         (i <= MIN(ulen,vlen + k))
	       ) 
	     )
      {
        while(i < ulen && j < vlen) 
        {
	  COMPARESYMBOLSSEP(i,j);
          i++;
          j++;
        }
      
        alwaysMININFINITYINT = false;
	tmpfront.dptabrow = i;
     
	// MAX(bigt_tmp, S(i+j,d));
        if (S (i + j, d) > bigt_tmp)
        {
          xdropbest->score = bigt_tmp = S (i + j, d);
          xdropbest->ivalue = (unsigned int) i;
          xdropbest->jvalue = (unsigned int) j;
        }

        /*
	DEBUG5 (3,
                "nach phase %ld:   d = %ld   k = %ld  (i,j) = (%ld,%ld)\n",
                (Showsint) d, (Showsint) d, (Showsint) k, (Showsint) i,
                (Showsint) j);
        */
      }
       
      else 
      {
        alwaysMININFINITYINT = false;
        tmpfront.dptabrow = 
	  fronts->spaceMyfrontvalue[ACCESSTOFRONT(d - 1, k)].dptabrow;
        /*
        DEBUG5 (3,
                "nach phase %ld:   d = %ld   k = %ld   (i,j) = (%ld,%ld)\n",
                (Showsint) d, (Showsint) d ,(Showsint) k, 
		(Showsint) tmpfront.dptabrow, (Showsint) tmpfront.dptabrow-k);
        */
      }

 
      STOREINARRAYFRONTS (fronts, ACCESSTOFRONT (d, k), Myfrontvalue,
                            tmpfront);

      /*
      DEBUG3(3, "fronts(%ld,%ld) = rowindex = %ld\n", 
	  (Showsint) d,
	  (Showsint) k,
	  (Showsint) fronts->spaceMyfrontvalue[ACCESSTOFRONT(d,k)].dptabrow);
      */
      
      // delete value for test for 
      // INSERTIONBIT/REPLACEMENTBIT/DELETIONBIT above
      tmpfront.dptabdirection = (unsigned char) 0;


    }   // end for loop

    // if all front values are MINUSINFINITYINT, 
    //   aligment prematurely finished if allowedMININFINITYINTgenerations 
    //   exceeded
    //  (full front has already ended at d - currentMININFINITYINTgeneration).
    if (alwaysMININFINITYINT)
    {
      currentMININFINITYINTgeneration++;
      if(currentMININFINITYINTgeneration > allowedMININFINITYINTgenerations)
      {
	d = d - currentMININFINITYINTgeneration;
        /*
	DEBUG1(1,
	    "Alignment dropped prematurely because of threshold "
	    "at distance = "
	    "%ld\n",
	    (Showsint) d);
        */
	break;
      }
    }
    else
    {
      currentMININFINITYINTgeneration = 0;
      alwaysMININFINITYINT = true;
    }
    
    /*
    DEBUG4 (3,
            "nach phase %ld:   best score so far = T' = %ld, (i,j) = "
	    "(%lu,%lu)\n\n",
            (Showsint) d, (Showsint) bigt_tmp,
            (Showuint) xdropbest->ivalue, (Showuint) xdropbest->jvalue);
    */
    STOREINARRAY (&big_t, Xdropscore, 10, bigt_tmp);

    // fill out of bounds values of MINUSINFINITYINT
    // needed for showmatrix function
    for (k = -d; k < lbound - 1; k++)
    {
      tmpfront.dptabrow = MINUSINFINITYINT;
      STOREINARRAYFRONTS (fronts, ACCESSTOFRONT (d, k), Myfrontvalue,
                          tmpfront);
    }
    for (k = ubound + 2; k <= d; k++)
    {
      tmpfront.dptabrow = MINUSINFINITYINT;
      STOREINARRAYFRONTS (fronts, ACCESSTOFRONT (d, k), Myfrontvalue,
                          tmpfront);
    }

    // alignment finished 
    if ((-d <= end_k && end_k <= d) &&
        fronts->spaceMyfrontvalue[ACCESSTOFRONT (d, end_k)].dptabrow == ulen)
    {
      /*
      DEBUG3 (2,
              "nach phase %ld: endpoint (%ld,%ld) reached. "
	      "alignment finished.\n",
              (Showsint) d, (Showsint) ulen, (Showsint) vlen);
      */
      break;
    }

    // pruning lower bound
    // lbound may decrease by one or increase/stays the same
    for (k = lbound - 1; k <= ubound + 1; k++)
    {
      if (fronts->spaceMyfrontvalue[ACCESSTOFRONT (d, k)].dptabrow >
          MINUSINFINITYINT)
      {
        lbound = k;
        break;
      }
    }

    //pruning upper bound
    // ubound may increase by one or decrease/stays the same
    for (k = ubound + 1; k >= lbound - 1; k--)
    {
      if (fronts->spaceMyfrontvalue[ACCESSTOFRONT (d, k)].dptabrow >
          MINUSINFINITYINT)
      {
        ubound = k;
        break;
      }
    }

    //handling boundaries lower bound
    lboundtmp = lbound;
    for (k = 0; k >= lbound; k--)
    {
      if (fronts->spaceMyfrontvalue[ACCESSTOFRONT (d, k)].dptabrow == vlen + k)
      {
        // debugging 12.11.06, vorher:
        //lboundtmp = k + 2; //bei Zhang so
        lboundtmp = k;
        break;
      }
    }

    //handling boundaries upper bound
    uboundtmp = ubound;
    for (k = 0; k <= ubound; k++)
    {
      if (fronts->spaceMyfrontvalue[ACCESSTOFRONT (d, k)].dptabrow == ulen)
      {
        // debugging 12.11.06, vorher:
        //uboundtmp = k - 2; //bei Zhang so
        uboundtmp = k;
        break;
      }
    }

    lbound = MAX (lbound, lboundtmp);
    ubound = MIN (ubound, uboundtmp);
  }

  /*
  DEBUG0(1, "Alignment results:\n");
  DEBUG1(1, "distance = %ld\n", (Showsint) d);
  DEBUG2(1,
           "best (prefix-) alignment score from (i,j) = (0,0) to "
	   "(%lu,%lu):\n",
           (Showuint) xdropbest->ivalue, (Showuint) xdropbest->jvalue);
  DEBUG1(1, "score T' = %ld\n", (Showsint) bigt_tmp);
  */

  FREEARRAY (&big_t, Xdropscore);
}

//#include "myxdrop.gen"
