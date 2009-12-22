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
#include "seqpos-def.h"
#include "verbose-def.h"

/*
  The following string is used to trigger the usage of gap costs
  for global chaining.
*/

#define GAPCOSTSWITCH        "gc"

/*
  The following string is used to trigger the use of a chaining algorithm
  allowing for overlaps between the hits.
*/

#define OVERLAPSWITCH        "ov"

/* the followin type is used position describing fragments
   to be chained */

typedef Seqpos GtChainpostype;

/* the following type is used for scores of chains */

typedef long GtChainscoretype;

/* the anonymous type for a chain */

typedef struct GtChain GtChain;

/* the anonymous type for a table storing the fragments */

typedef struct GtFragmentinfotable GtFragmentinfotable;

/* the following type is used for relevant values to output for a chained
   frgament */

typedef struct
{
  GtChainpostype startpos[2], /* start of fragments in the 2 dimensions,
                                 userdef */
                 endpos[2];  /* end of fragments in the 2 dimensions, userdef */
  GtChainscoretype weight; /* weight of fragment, user defined */
} GtFragmentvalues;

/*
  the type of function to report chains.
*/

typedef void (*GtChainprocessor)(void *,
                                 const GtFragmentinfotable *,
                                 const GtChain *);

/* the type of value describing how to chain */

typedef struct GtChainmode GtChainmode;

/* the constructor for tables of fragments */

GtFragmentinfotable *gt_chain_fragmentinfotable_new(
                           unsigned long numberoffragments);

/* the destructor for tables of fragments */

void gt_chain_fragmentinfotable_delete(GtFragmentinfotable *fragmentinfotable);

/* the function for emptying tables of fragments, without freeing the space */

void gt_chain_fragmentinfotable_empty(GtFragmentinfotable *fragmentinfotable);

/* the following function adds the relevant values describing a fragment */

void gt_chain_fragmentinfotable_add(GtFragmentinfotable *fragmentinfotable,
                                    const GtFragmentvalues *infragment);

/* the following functions reads a file describing fragments in open format.
   It returns the corresponding fragment table */

GtFragmentinfotable *gt_chain_analyzeopenformatfile(double weightfactor,
                                                    const char *matchfile,
                                                    GtError *err);

/* the function to fill the gap values for all fragments */

void gt_chain_fillthegapvalues(GtFragmentinfotable *fragmentinfotable);

/* the function to sort an array of fragments */

void gt_chain_possiblysortfragments(Verboseinfo *verboseinfo,
                                    GtFragmentinfotable *fragmentinfotable,
                                    unsigned int presortdim);

/* the constructor for chainmode objects. Use err = NULL to print
   error messages to stderr. */

GtChainmode *gt_chain_chainmode_new(unsigned long maxgap,
                                    bool globalset,
                                    const char *globalargs,
                                    bool localset,
                                    const char *localargs,
                                    GtError *err);

/* the destructor for chainmode objects */

void gt_chain_chainmode_delete(GtChainmode *gtchainmode);

/* the constructor for chains */

GtChain *gt_chain_chain_new(void);

/* the descructor for chains */

void gt_chain_chain_delete(GtChain *chain);

/* the function to perform the fast chaining algorithms */

void gt_chain_fastchaining(const GtChainmode *chainmode,
                           GtChain *chain,
                           GtFragmentinfotable *fragmentinfotable,
                           bool gapsL1,
                           unsigned int presortdim,
                           bool withequivclasses,
                           GtChainprocessor chainprocessor,
                           void *cpinfo,
                           Verboseinfo *verboseinfo);

/* obtain the score of a chain */

GtChainscoretype gt_chain_chainscore(const GtChain *chain);

/* obtain the length of a chain */

unsigned long gt_chain_chainlength(const GtChain *chain);

/* store the values of element idx in given chain in the first parameter */

void gt_chain_extractchainelem(GtFragmentvalues *value,
                               const GtFragmentinfotable *fragmentinfotable,
                               const GtChain *chain,
                               unsigned long idx);

/* print a chain element to the given file pointer */

void gt_chain_printchainelem(FILE *outfp,const GtFragmentvalues *value);

#endif
