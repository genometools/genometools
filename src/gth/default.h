/*
  Copyright (c) 2003-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef DEFAULT_H
#define DEFAULT_H

#include "core/translator_api.h"

/*
  This file contains all default values of GenomeThreader.
  Default values can be changed by command line options.
*/

/* default score matrix */
#define DEFAULT_SCOREMATRIX                     "BLOSUM62"

/* defaults for data preprocessing options */
#define DEFAULT_PROTEINSMAP                     "protein"

/* default translation table */
#define DEFAULT_TRANSLATIONTABLE                GT_TRANSLATOR_STANDARD_SCHEME

/* default minimum match length for (new) similarity filter */
#define DEFAULT_MINMATCHLENGTH                  20
#define DEFAULT_SEEDLENGTH                      18
#define DEFAULT_EXDROP                          2

/* default vmatch parameter for protein matching (taken from GS2) */
#define DEFAULT_PRMINMATCHLEN                   24
#define DEFAULT_PRSEEDLENGTH                    10
#define DEFAULT_PRHDIST                         4

#define DEFAULT_ONLINE                          false
#define DEFAULT_INVERSE                         false
#define DEFAULT_EXACT                           false
#define DEFAULT_EDIST                           false
#define DEFAULT_NOAUTOINDEX                     false
#define DEFAULT_CREATEINDICESONLY               false
#define DEFAULT_SKIPINDEXCHECK                  false
#define DEFAULT_MASKPOLYATAILS                  false
#define DEFAULT_MAXNUMOFMATCHES                 0

#define DEFAULT_FRAGWEIGHTFACTOR                0.5
#define DEFAULT_GCMAXGAPWIDTH                   1000000
#define DEFAULT_RARE                            0
#define DEFAULT_PARALOGS                        false
#define DEFAULT_ENRICHCHAINS                    false
#define DEFAULT_GCMINCOVERAGE                   50
#define DEFAULT_STOPAFTERCHAINING               false

/* default values for the (generic) splice site model */
#define DEFAULT_DISABLEU12INTRONMODEL           false

#define DEFAULT_GENERIC_GT_DONORPROB            ((GthFlt) 0.05)
#define DEFAULT_NONGENERIC_GT_DONORPROB         ((GthFlt) 0.00005)
#define DEFAULT_GENERIC_GC_DONORPROB            ((GthFlt) 0.002)
#define DEFAULT_NONGENERIC_GC_DONORPROB         ((GthFlt) 0.00002)
#define DEFAULT_GENERIC_AT_DONORPROB            ((GthFlt) 0.002)
#define DEFAULT_NONGENERIC_AT_DONORPROB         ((GthFlt) 0.00002)

#define DEFAULT_GENERIC_AG_ACCEPTORPROB         ((GthFlt) 0.05)
#define DEFAULT_NONGENERIC_AG_ACCEPTORPROB      ((GthFlt) 0.00005)
#define DEFAULT_GENERIC_AC_ACCEPTORPROB         ((GthFlt) 0.002)
#define DEFAULT_NONGENERIC_AC_ACCEPTORPROB      ((GthFlt) 0.00002)

#define DEFAULT_GENERIC_OTHERSPLICESITEPROB     ((GthFlt) 0.0001)
#define DEFAULT_NONGENERIC_OTHERSPLICESITEPROB  ((GthFlt) 0.000001)

#define DEFAULT_U12_TYPEDONORPROB               ((GthFlt) 0.99)
#define DEFAULT_U12_TYPEDONORPROBONEMISMATCH    ((GthFlt) 0.9)

/* default values for intron cutouts */
#define DEFAULT_INTRONCUTOUT         false
#define DEFAULT_AUTOICMAXMATRIXSIZE  0
#define DEFAULT_MAXMATRIXSIZE        UNDEF_ULONG
#define DEFAULT_ICINITIALDELTA       50
#define DEFAULT_ICITERATIONS         2
#define DEFAULT_ICDELTAINCREASE      50
#define DEFAULT_ICMINREMLENGTH       10

#define GTH_DEFAULT_NOICININTRONCHECK    false
#define GTH_DEFAULT_FREEINTRONTRANS      false
#define GTH_DEFAULT_DPMINEXONLENGTH      5
#define GTH_DEFAULT_DPMININTRONLENGTH    50
#define GTH_DEFAULT_SHORTEXONPENALTY     100.0
#define GTH_DEFAULT_SHORTINTRONPENALTY   100.0

#define GTH_DEFAULT_PROBIES              0.5
#define GTH_DEFAULT_PROBDELGEN           0.03
#define GTH_DEFAULT_IDENTITYWEIGHT       2.0
#define GTH_DEFAULT_MISMATCHWEIGHT       -2.0
#define GTH_DEFAULT_UNDETCHARWEIGHT      0.0
#define GTH_DEFAULT_DELETIONWEIGHT       -5.0
#define GTH_DEFAULT_WZEROTRANSITION      80
#define GTH_DEFAULT_WDECREASEDOUTPUT     80
#define GTH_DEFAULT_DETECTSMALLEXONS     false

#define GTH_DEFAULT_CUTOFFSMINEXONLEN    5
#define GTH_DEFAULT_SCOREMINEXONLEN      50

#define DEFAULT_MINAVERAGESSP        0.5       /* minimum average splice site
                                                  probability */
#define DEFAULT_MIN_ALIGNMENTSCORE   0.0       /* minimum aligmment score */
#define DEFAULT_MAX_ALIGNMENTSCORE   1.0       /* maximum aligmment score */
#define DEFAULT_MIN_COVERAGE         0.0       /* minimum coverage */

/* the coverage computation can produce a coverage larger than 1.0
   and Volker Brendel thinks it is best to leave it this way */
#define DEFAULT_MAX_COVERAGE         9999.99   /* maximum coverage */

#define DEFAULT_INTERMEDIATE         false     /* stop after SA calculation */

/* defaults for sorting of AGSs */
#define DEFAULT_SORTAGS              false
#define DEFAULT_SORTAGSWF            1.0

#define DEFAULT_DISABLECLUSTERSAS    true      /* XXX !!! */

#define DEFAULT_FIRSTALSHOWN         0         /* default for firstalshown */
                                               /* default for showfullintrons */
#define DEFAULT_SHOWINTRONMAXLEN DOUBLEALIGNMENTLINEWIDTH

#define DEFAULT_VERBOSESEQS          false
#define DEFAULT_SKIPALIGNMENTOUT     false
#define DEFAULT_SHOWSEQNUMS          false
#define DEFAULT_SHOWMATUREMRNAS      false
#define DEFAULT_XMLOUT               false

#define DEFAULT_MINORFLENGTH         64         /* default for minORFlength */

#define DEFAULT_GS2OUT               false      /* default for gs2out */

#define DEFAULT_EXONDISTRI           false
#define DEFAULT_INTRONDISTRI         false
#define DEFAULT_MATCHNUMDISTRI       false
#define DEFAULT_REFSEQCOVDISTRI      false

#define DEFAULT_COMMENTS             false
#define DEFAULT_CHAINWITHOVERLAPS    false
#define DEFAULT_PROTEINEXONPENAL     false
#define DEFAULT_PROTEINDPUPGRADE     false

/* default for test options */
#define DEFAULT_SHOWEOPS             false

#endif
