/*
  Copyright (c) 2003-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2007 Center for Bioinformatics, University of Hamburg

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

#ifndef GTHSIMFILTERPARAM_H
#define GTHSIMFILTERPARAM_H

#include "gth/gthintroncutout.h"

typedef struct {
  /* dna vmatch call: */
  unsigned long minmatchlength,  /* minimum match length */
                seedlength,      /* seed length */
                exdrop,          /* xdrop value for edit distance extension */

  /* protein vmatch call: */
                prminmatchlen,   /* minimum match length */
                prseedlength,    /* seed length */
                prhdist,         /* haming distance */

  /* other stuff: */
                maxnumofmatches, /* the maximum number of matches (per genomic
                                    file, per reference sequence) */
                rare;            /* maximum number of (rare) matches */
  bool online,                   /* use online algorithm */
       inverse,                  /* invert query and index */
       exact,                    /* compute exact matches */
       edist,                    /* use edist instead of exdrop */
       noautoindex,              /* do not create indices automatically */
       createindicesonly,        /* stop the program flow after the indices have
                                    been created */
       skipindexcheck,           /* skip index check (in preprocessing phase) */
       maskpolyAtails,           /* create and use masked files for vmatch call
                                  */
       paralogs,                 /* compute paralogous genes
                                    (different chaining procedure) */
       enrichchains,             /* enrich chains with additional matches */
       stopafterchaining;        /* stop gth after chaining phase */
  bool jump_table;               /* use jump table in DP */
  Introncutoutinfo introncutoutinfo; /* parameter for intron cutout */
} Gthsimfilterparam;

#endif
