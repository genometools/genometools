/* Copyright (C) 2007 by David Ellinghaus <dellinghaus@zbh.uni-hamburg.de> */

#ifndef REPEATTYPES_H
#define REPEATTYPES_H

#include <stdbool.h>

//#include "types.h"
#include "libgtmatch/arraydef.h"
#include "libgtmatch/seqpos-def.h"
#include "libgtmatch/sarr-def.h"
//#include "virtualdef.h" ersetzen

/*
 The datatype Repeat stores information about the maximal repeats (seeds).
 */
typedef struct
{
Seqpos pos1;         // first position of maximal repeat (seed)
Seqpos offset;       // second position = pos1 + offset
Seqpos len;          // length of maximal repeat
unsigned long contignumber; // number of contig for this repeat
} Repeat;   // \Typedef{Repeat}

DECLAREARRAYSTRUCT(Repeat);

/*
 The datatype RepeatInfo stores all maximal repeats (seeds) and
 information about the length and distance constraints.
 */
typedef struct
{
ArrayRepeat repeats; // array of maximal repeats (seeds)
unsigned long lmin;        // minimum allowed length of a LTR
unsigned long lmax;        // maximum allowed length of a LTR
unsigned long dmin;        // minimum distance between LTRs
unsigned long dmax;        // maximum distance between LTRs
Suffixarray *suffixarrayptr; 
/* pointer on virtualtree, is needed in function simpleexactselfmatchstore, repeats.c */ 
//ersetzen!!! 
} RepeatInfo;   // \Typedef{RepeatInfo}

/*
 The datatype SubRepeatInfo stores information about the maximal repeats
 for the TSD detection.
 */
typedef struct
{
ArrayRepeat repeats; // array of maximal repeats for TSDs
unsigned longlmin;           // minimal length of TSD
unsigned longlmax;           // maximal length of TSD
Seqpos separatorpos;   // position of separator of the 
                     // two concatenated sequences 
Seqpos offset1;        // offset1 for absolute position 1 in sequence
Seqpos offset2;        // offset2 for absolute position 2 in sequence
                     // pos1 < pos2
} SubRepeatInfo; // \Typedef{SubRepeatInfo}

/*
 The datatype LTRboundaries stores all information of one predicted
 LTR element.
 */
typedef struct
{
  unsigned long contignumber; // number of sequence in multiseq of virtualtree
  Seqpos leftLTR_5,    // 5' boundary of left LTR
         leftLTR_3,    // 3' boundary of left LTR
         rightLTR_5,   // 5' boundary of right LTR 
         rightLTR_3;   // 3' boundary of right LTR
  Seqpos lenleftTSD,
         lenrightTSD;
  bool tsd,            // If true, then TSDs exist.
       motif_near_tsd, // If true, then motif near the TSD exists.
       motif_far_tsd,  /* If true, then motif at the inner borders of
                          LTRs exist. */
       lengthdistconstraint; /* If true, length and distance constraints
                                are fulfilled */
  double similarity;   // similarity value of LTRs
  bool skipped;        /* if skipped then because of an overlap
                          with a higher similarity prediction or
                          because of "noclusterallowed" option */
} LTRboundaries;                 // \Typedef{LTRboundaries}

DECLAREARRAYSTRUCT(LTRboundaries);
#endif
