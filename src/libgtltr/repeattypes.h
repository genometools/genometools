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

#ifndef REPEATTYPES_H
#define REPEATTYPES_H

#include <stdbool.h>

#include "libgtcore/arraydef.h"
#include "libgtmatch/esa-seqread.h"
#include "libgtmatch/sarr-def.h"

/* The datatype Repeat stores information about the maximal repeats (seeds).*/
typedef struct
{
  Seqpos pos1;         /* first position of maximal repeat (seed) */
  Seqpos offset;       /* second position = pos1 + offset */
  Seqpos len;          /* length of maximal repeat  */
  unsigned long contignumber; /* number of contig for this repeat */
} Repeat;

DECLAREARRAYSTRUCT(Repeat);

/* The datatype RepeatInfo stores all maximal repeats (seeds) and */
/* information about the length and distance constraints. */
typedef struct
{
  ArrayRepeat repeats; /* array of maximal repeats (seeds) */
  unsigned long lmin;        /* minimum allowed length of a LTR */
  unsigned long lmax;        /* maximum allowed length of a LTR */
  unsigned long dmin;        /* minimum distance between LTRs */
  unsigned long dmax;        /* maximum distance between LTRs */
  Sequentialsuffixarrayreader *ssarptr;
/* pointer on suffixarray, is needed in function simpleexactselfmatchstore,
   repeats.c */
} RepeatInfo;

/* The datatype SubRepeatInfo stores information about the maximal repeats */
/* for the TSD detection. */
typedef struct
{
  ArrayRepeat repeats; /* array of maximal repeats for TSDs */
  unsigned long lmin;   /* minimal length of TSD */
  unsigned long lmax;   /* maximal length of TSD */
  Seqpos offset1;      /* offset1 for absolute position 1 in sequence */
  Seqpos offset2;      /* offset2 for absolute position 2 in sequence */
                       /* pos1 < pos2 */
} SubRepeatInfo;

/* The datatype LTRboundaries stores all information of one predicted */
/* LTR element. */
typedef struct
{
  unsigned long contignumber; /* number of sequence in multiseq of
                                 virtualtree */
  Seqpos leftLTR_5,    /* 5' boundary of left LTR */
         leftLTR_3,    /* 3' boundary of left LTR */
         rightLTR_5,   /* 5' boundary of right LTR */
         rightLTR_3;   /* 3' boundary of right LTR */
  Seqpos lenleftTSD,
         lenrightTSD;
  bool tsd,            /* If true, then TSDs exist. */
       motif_near_tsd, /* If true, then motif near the TSD exists. */
       motif_far_tsd,  /* If true, then motif at the inner borders of  */
                       /* LTRs exist. */
       lengthdistconstraint; /* If true, length and distance constraints */
                             /* are fulfilled */
  double similarity;   /* similarity value of LTRs */
  bool skipped;        /* if skipped then because of an overlap
                          with a higher similarity prediction or
                          because of "noclusterallowed" option */
} LTRboundaries;

DECLAREARRAYSTRUCT(LTRboundaries);

/* The datatype Motif stores information about the specified motif. */
typedef struct
{
  Str *str_motif;
  Uchar firstleft, /* first character of left motif instance */
        secondleft,    /* second character of left motif instance */
        firstright,    /* first character of right motif instance */
        secondright;   /* second character of right motif instance */
  unsigned int allowedmismatches; /* number of allowed mismatches in the four */
                                  /*character motif */
} Motif;
#endif
