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

#ifndef CHAINDEF_H
#define CHAINDEF_H

#include "core/arraydef.h"
#include "core/error_api.h"
#include "core/str_api.h"
#include "seqpos-def.h"
#include "verbose-def.h"

typedef long GtChainscoretype;

/*
  The following type defines the possible kinds of chaining.
  The mode can be one of the two following values.
*/

typedef enum
{
  GLOBALCHAINING,              /* global chaining without gap costs */
  GLOBALCHAININGWITHGAPCOST,   /* global chaining with L1 gap costs */
  GLOBALCHAININGWITHOVERLAPS,  /* chaining with overlaps */
  LOCALCHAININGMAX,            /* local chaining; one maximum is reported */
  LOCALCHAININGTHRESHOLD,      /* local chaining; all chains >= minscore */
  LOCALCHAININGBEST,           /* local chaining; k best local chains */
  LOCALCHAININGPERCENTAWAY     /* local chaining; percent away from best */
} GtChainkind;

/*
  A chain consists of an array of integers. These refer to the array of
  fragment informations.
*/

typedef unsigned long GtChainref;

GT_DECLAREARRAYSTRUCT(GtChainref);

typedef struct
{
  GtArrayGtChainref chainedfragments;
  GtChainscoretype scoreofchain;
} GtChain;

/*
  We use functions of the following type to report chains.
*/

typedef int (*GtChainprocessor)(void *,GtChain *,GtError *err);

/*
  The following type defines the chain mode consisting of a chainkind.
  If chainkind = LOCALCHAININGTHRESHOLD, then an additional
  component minimumscore is used.
  If chainkind = LOCALCHAININGBEST, then  an additional
  component howmanybest is used.
  If chainkind = LOCALCHAININGPERCENTAWAY, then  an additional
  component percentawayfrombest is defined
*/

typedef Seqpos GtChainpostype;

typedef struct
{
  GtChainkind chainkind;
  GtChainpostype maxgapwidth;  /* 0 if undefined or
                                  otherwise maximal width of gap */
  GtChainscoretype minimumscore; /* only defined if
                                  chainkind = LOCALCHAININGTHRESHOLD */
  unsigned long howmanybest,   /* only defined if
                                  chainkind = LOCALCHAININGBEST */
                percentawayfrombest;  /* only defined if
                                         chainkind = LOCALCHAININGPERCENTAWAY */
} GtChainmode;

typedef struct
{
  bool silent;
  GtChainmode chainmode;
  double weightfactor;
  const GtStr *matchfile,
              *outprefix;
  Verboseinfo *verboseinfo;
} Chaincallinfo;

typedef struct GtFragmentinfotable GtFragmentinfotable;

GtFragmentinfotable *gt_chain_fragmentinfotable_new(
                           unsigned long numberoffragments);

void gt_chain_fragmentinfotable_delete(GtFragmentinfotable *fragmentinfotable);

void gt_chain_fragmentinfotable_add(GtFragmentinfotable *fragmentinfotable,
                                    GtChainpostype start1,
                                    GtChainpostype end1,
                                    GtChainpostype start2,
                                    GtChainpostype end2,
                                    GtChainscoretype weight);

void gt_chain_fillthegapvalues(GtFragmentinfotable *fragmentinfotable);

int gt_chain_fastchaining(const GtChainmode *chainmode,
                    GtChain *chain,
                    GtFragmentinfotable *fragmentinfotable,
                    bool gapsL1,
                    unsigned int presortdim,
                    bool withequivclasses,
                    GtChainprocessor chainprocessor,
                    void *cpinfo,
                    Verboseinfo *verboseinfo,
                    GtError *err);

void gt_chain_possiblysortopenformatfragments(
                             Verboseinfo *verboseinfo,
                             GtFragmentinfotable *fragmentinfotable,
                             unsigned int presortdim);

#endif
