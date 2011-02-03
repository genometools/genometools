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

#ifndef CALL_INFO_H
#define CALL_INFO_H

#include "core/array.h"
#include "gth/gthoutput.h"
#include "gth/dp_options_core.h"
#include "gth/dp_options_est.h"
#include "gth/dp_options_postpro.h"
#include "gth/gthalphatype.h"
#include "gth/gthsimfilterparam.h"
#include "gth/sa_filter.h"
#include "gth/splice_site_model.h"

/* contains the all the parameters, which are passed to the program */
typedef struct {
  unsigned int translationtable,     /* translation table used at various places
                                      */
               speciesnum,           /* the species number determining the model
                                      */
               firstalshown,         /* number of cDNA/EST alignments shown
                                       (GS2=maxnest) */
               gcmaxgapwidth,        /* maximum gap width for global chains */
               gcmincoverage;        /* minimum coverage of global chains
                                        regarding to the reference sequence */
  char *progname;                    /* name of this binary (e.g., ``gth'') */
  GtStr *scorematrixfile;            /* file name of amino acid substitution
                                        matrix */
  Gthsimfilterparam simfilterparam;  /* the parameter for the similarity filter
                                      */
  GthDPOptionsCore *dp_options_core; /* the core DP options */
  GthDPOptionsEST *dp_options_est;   /* the cDNA/EST DP options */
  GthDPOptionsPostpro *dp_options_postpro ;/* DP post processing options for
                                              processing ``raw'' spliced
                                              alignments */
  double fragweightfactor,           /* fragment weight factor */
         minaveragessp;              /* minimum average splice site probability
                                      */
  bool intermediate,                 /* stop after SA computation */
       proteinexonpenal,             /* add short exon penalty in protein DP */
       disableclustersas,            /* disable the clustering of SAs in
                                        consensus phase */
       cdnaforwardonly;
  GthSAFilter *sa_filter;            /* the spliced alignment filter */
  GthDuplicateCheck duplicate_check; /* the modus use for duplicate checks */
  GthSpliceSiteModel *splice_site_model; /* the splice site model */
  GthOutput *out;                    /* bundles all output related variables */
} GthCallInfo;

GthCallInfo*  gth_call_info_new(const char *progname);
void          gth_call_info_delete(GthCallInfo*);

#endif
