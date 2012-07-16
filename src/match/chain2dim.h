/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#ifndef CHAIN2DIM_H
#define CHAIN2DIM_H

#include "core/error_api.h"
#include "core/logger.h"

/*
  The following string is used to trigger the usage of gap costs
  for global chaining.
*/

#define GT_CHAIN2DIM_GAPCOSTSWITCH        "gc"

/*
  The following string is used to trigger the use of a chaining algorithm
  allowing for overlaps between the hits.
*/

#define GT_CHAIN2DIM_OVERLAPSWITCH        "ov"

#define GT_CHAIN2DIM_ALLSWITCH            "all"

/* the followin type is used for the position values in the matches to be
   chained */

typedef unsigned long GtChain2Dimpostype;

/* the following type is used for scores of chains */

typedef long GtChain2Dimscoretype;

/* the anonymous type for a chain */

typedef struct GtChain2Dim GtChain2Dim;

/* the anonymous type for a table storing the matches */

typedef struct GtChain2Dimmatchtable GtChain2Dimmatchtable;

/* the following type is used for relevant values to output for a chained
   frgament */

typedef struct
{
  GtChain2Dimpostype startpos[2], /* start of matches in the 2 dimensions,
                                 userdef */
                 endpos[2];  /* end of matches in the 2 dimensions, userdef */
  GtChain2Dimscoretype weight; /* weight of match, user defined */
} GtChain2Dimmatchvalues;

/*
  the type of function to report chains.
*/

typedef void (*GtChain2Dimprocessor)(void *,
                                 const GtChain2Dimmatchtable *,
                                 const GtChain2Dim *);

/* the type of value describing how to chain */

typedef struct GtChain2Dimmode GtChain2Dimmode;

/* the constructor for tables of matches */

GtChain2Dimmatchtable *gt_chain_matchtable_new(unsigned long numberofmatches);

/* the destructor for tables of matches */

void gt_chain_matchtable_delete(GtChain2Dimmatchtable *matchtable);

/* the function for emptying a table of matches, without freeing the space */

void gt_chain_matchtable_empty(GtChain2Dimmatchtable *matchtable);

/* the following function adds the relevant values describing a match */

void gt_chain_matchtable_add(GtChain2Dimmatchtable *matchtable,
                                    const GtChain2Dimmatchvalues *inmatch);

/* the following functions reads a file describing matches in open format.
   It returns the corresponding table of matches. */

GtChain2Dimmatchtable *gt_chain_analyzeopenformatfile(double weightfactor,
                                                  const char *matchfile,
                                                  GtError *err);

/* the function to fill the gap values for all matches */

void gt_chain_fillthegapvalues(GtChain2Dimmatchtable *matchtable);

/* function to apply an additional weight to the elements to be chained */

void gt_chain_applyweight(double weightfactor,
                          GtChain2Dimmatchtable *matchtable);

/* the function to sort an array of matches */

void gt_chain_possiblysortmatches(GtLogger *logger,
                                  GtChain2Dimmatchtable *matchtable,
                                  unsigned int presortdim);

/* the constructor for chainmode objects. Use err = NULL to print
   error messages to stderr. */

GtChain2Dimmode *gt_chain_chainmode_new(unsigned long maxgap,
                                        bool globalset,
                                        const char *globalargs,
                                        bool localset,
                                        const char *localargs,
                                        GtError *err);

/* the destructor for chainmode objects */

void gt_chain_chainmode_delete(GtChain2Dimmode *chainmode);

/* the constructor for chains */

GtChain2Dim *gt_chain_chain_new(void);

/* the descructor for chains */

void gt_chain_chain_delete(GtChain2Dim *chain);

/* the function to perform the fast chaining algorithms */

void gt_chain_fastchaining(const GtChain2Dimmode *chainmode,
                           GtChain2Dim *chain,
                           GtChain2Dimmatchtable *matchtable,
                           bool gapsL1,
                           unsigned int presortdim,
                           bool withequivclasses,
                           GtChain2Dimprocessor chainprocessor,
                           void *cpinfo,
                           GtLogger *logger);

/* obtain the score of a chain */

GtChain2Dimscoretype gt_chain_chainscore(const GtChain2Dim *chain);

/* obtain the length of a chain */

unsigned long gt_chain_chainlength(const GtChain2Dim *chain);

/* return true iff chain is stored in reverse order */

bool gt_chain_storedinreverseorder(const GtChain2Dim *chain);

/* store the values of element idx in given chain in the first parameter */

void gt_chain_extractchainelem(GtChain2Dimmatchvalues *value,
                               const GtChain2Dimmatchtable *matchtable,
                               const GtChain2Dim *chain,
                               unsigned long idx);

/* print a chain element to the given file pointer */

void gt_chain_printchainelem(FILE *outfp,const GtChain2Dimmatchvalues *value);

#endif
