/*
  Copyright (c) 2003-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/trans_table_api.h"

/*
  This file contains all default values of GenomeThreader.
  Default values can be changed by command line options.
*/

/* default score matrix */
#define GTH_DEFAULT_SCOREMATRIX        "BLOSUM62"

/* defaults for data preprocessing options */
#define GTH_DEFAULT_PROTEINSMAP        "protein"

/* default translation table */
#define GTH_DEFAULT_TRANSLATIONTABLE   GT_STANDARD_TRANSLATION_SCHEME

/* default minimum match length for (new) similarity filter */
#define GTH_DEFAULT_MINMATCHLENGTH     20
#define GTH_DEFAULT_SEEDLENGTH         18
#define GTH_DEFAULT_EXDROP             2

/* default vmatch parameter for protein matching (taken from GS2) */
#define GTH_DEFAULT_PRMINMATCHLEN      24
#define GTH_DEFAULT_PRSEEDLENGTH       10
#define GTH_DEFAULT_PRHDIST            4

#define GTH_DEFAULT_ONLINE             false
#define GTH_DEFAULT_INVERSE            false
#define GTH_DEFAULT_EXACT              false
#define GTH_DEFAULT_EDIST              false
#define GTH_DEFAULT_NOAUTOINDEX        false
#define GTH_DEFAULT_CREATEINDICESONLY  false
#define GTH_DEFAULT_SKIPINDEXCHECK     false
#define GTH_DEFAULT_MASKPOLYATAILS     false
#define GTH_DEFAULT_MAXNUMOFMATCHES    0

#define GTH_DEFAULT_FRAGWEIGHTFACTOR   0.5
#define GTH_DEFAULT_GCMAXGAPWIDTH      1000000
#define GTH_DEFAULT_RARE               0
#define GTH_DEFAULT_PARALOGS           false
#define GTH_DEFAULT_ENRICHCHAINS       false
#define GTH_DEFAULT_GCMINCOVERAGE      50
#define GTH_DEFAULT_STOPAFTERCHAINING  false

/* default values for the (generic) splice site model */
#define GTH_DEFAULT_DISABLEU12INTRONMODEL           false

#define GTH_DEFAULT_GENERIC_GT_DONORPROB            ((GthFlt) 0.05)
#define GTH_DEFAULT_NONGENERIC_GT_DONORPROB         ((GthFlt) 0.00005)
#define GTH_DEFAULT_GENERIC_GC_DONORPROB            ((GthFlt) 0.002)
#define GTH_DEFAULT_NONGENERIC_GC_DONORPROB         ((GthFlt) 0.00002)
#define GTH_DEFAULT_GENERIC_AT_DONORPROB            ((GthFlt) 0.002)
#define GTH_DEFAULT_NONGENERIC_AT_DONORPROB         ((GthFlt) 0.00002)

#define GTH_DEFAULT_GENERIC_AG_ACCEPTORPROB         ((GthFlt) 0.05)
#define GTH_DEFAULT_NONGENERIC_AG_ACCEPTORPROB      ((GthFlt) 0.00005)
#define GTH_DEFAULT_GENERIC_AC_ACCEPTORPROB         ((GthFlt) 0.002)
#define GTH_DEFAULT_NONGENERIC_AC_ACCEPTORPROB      ((GthFlt) 0.00002)

#define GTH_DEFAULT_GENERIC_OTHERSPLICESITEPROB     ((GthFlt) 0.0001)
#define GTH_DEFAULT_NONGENERIC_OTHERSPLICESITEPROB  ((GthFlt) 0.000001)

#define GTH_DEFAULT_U12_TYPEDONORPROB               ((GthFlt) 0.99)
#define GTH_DEFAULT_U12_TYPEDONORPROBONEMISMATCH    ((GthFlt) 0.9)

#define GTH_DEFAULT_JUMPTABLE            false

/* default values for intron cutouts */
#define GTH_DEFAULT_INTRONCUTOUT         false
#define GTH_DEFAULT_AUTOICMAXMATRIXSIZE  0
#define GTH_DEFAULT_MAXMATRIXSIZE        UNDEF_ULONG
#define GTH_DEFAULT_ICINITIALDELTA       50
#define GTH_DEFAULT_ICITERATIONS         2
#define GTH_DEFAULT_ICDELTAINCREASE      50
#define GTH_DEFAULT_ICMINREMLENGTH       10

#define GTH_DEFAULT_NOICININTRONCHECK    false
#define GTH_DEFAULT_FREEINTRONTRANS      false
#define GTH_DEFAULT_DPMINEXONLENGTH      5
#define GTH_DEFAULT_DPMININTRONLENGTH    50
#define GTH_DEFAULT_SHORTEXONPENALTY     100.0
#define GTH_DEFAULT_SHORTINTRONPENALTY   100.0

#define GTH_DEFAULT_JTOVERLAP            5
#define GTH_DEFAULT_JTDEBUG              false

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

#define GTH_DEFAULT_MINAVERAGESSP        0.5      /* minimum average splice site
                                                     probability */
#define GTH_DEFAULT_MIN_ALIGNMENTSCORE   0.0      /* minimum aligmment score */
#define GTH_DEFAULT_MAX_ALIGNMENTSCORE   1.0      /* maximum aligmment score */
#define GTH_DEFAULT_MIN_COVERAGE         0.0      /* minimum coverage */

/* the coverage computation can produce a coverage larger than 1.0
   and Volker Brendel thinks it is best to leave it this way */
#define GTH_DEFAULT_MAX_COVERAGE         9999.99  /* maximum coverage */

#define GTH_DEFAULT_INTERMEDIATE         false    /* stop after SA
                                                     calculation */

/* defaults for sorting of AGSs */
#define GTH_DEFAULT_SORTAGS              false
#define GTH_DEFAULT_SORTAGSWF            1.0

#define GTH_DEFAULT_DISABLECLUSTERSAS    true     /* XXX !!! */

#define GTH_DEFAULT_FIRSTALSHOWN         0        /* default for firstalshown */

/* default for showfullintrons */
#define GTH_DEFAULT_SHOWINTRONMAXLEN     DOUBLEALIGNMENTLINEWIDTH

#define GTH_DEFAULT_VERBOSESEQS          false
#define GTH_DEFAULT_SKIPALIGNMENTOUT     false
#define GTH_DEFAULT_SHOWSEQNUMS          false
#define GTH_DEFAULT_SHOWMATUREMRNAS      false
#define GTH_DEFAULT_XMLOUT               false

#define GTH_DEFAULT_MINORFLENGTH         64       /* default for minORFlength */
#define GTH_DEFAULT_START_CODON          false
#define GTH_DEFAULT_FINAL_STOP_CODON     false

#define GTH_DEFAULT_GS2OUT               false    /* default for gs2out */
#define GTH_DEFAULT_MD5IDS               false

#define GTH_DEFAULT_EXONDISTRI           false
#define GTH_DEFAULT_INTRONDISTRI         false
#define GTH_DEFAULT_MATCHNUMDISTRI       false
#define GTH_DEFAULT_REFSEQCOVDISTRI      false

#define GTH_DEFAULT_COMMENTS             false
#define GTH_DEFAULT_CHAINWITHOVERLAPS    false
#define GTH_DEFAULT_PROTEINEXONPENAL     false
#define GTH_DEFAULT_PROTEINDPUPGRADE     false

/* default for test options */
#define GTH_DEFAULT_SHOWEOPS             false

#endif
