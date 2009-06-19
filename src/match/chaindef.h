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
#include "seqpos-def.h"

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

typedef int (*Chainprocessor)(void *,GtChain *);

/*
  The following type defines the chain mode consisting of a chainkind.
  If chainkind = LOCALCHAININGTHRESHOLD, then an additional
  component minimumscore is used.
  If chainkind = LOCALCHAININGBEST, then  an additional
  component howmanybest is used.
  If chainkind = LOCALCHAININGPERCENTAWAY, then  an additional
  component percentawayfrombest is defined
*/

typedef struct
{
  GtChainkind chainkind;
  Seqpos maxgapwidth;   /* 0 if undefined or otherwise maximal width of gap */
  GtChainscoretype minimumscore; /* only defined if
                                  chainkind = LOCALCHAININGTHRESHOLD */
  unsigned long howmanybest,   /* only defined if
                                  chainkind = LOCALCHAININGBEST */
                percentawayfrombest;  /* only defined if
                                         chainkind = LOCALCHAININGPERCENTAWAY */
} Chainmode;

typedef struct GtFragmentinfotable GtFragmentinfotable;

#endif
