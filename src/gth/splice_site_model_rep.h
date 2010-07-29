/*
  Copyright (c) 2005-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef SPLICE_SITE_MODEL_REP_H
#define SPLICE_SITE_MODEL_REP_H

#include "gth/splice_site_model.h"

struct GthSpliceSiteModel {
  /* the probabilities */
  GthFlt genericGTdonorprob,          /*     generic prob. of GT donor */
         nongenericGTdonorprob,       /* non-generic prob. of GT donor */
         genericGCdonorprob,          /*     generic prob. of GC donor */
         nongenericGCdonorprob,       /* non-generic prob. of GC donor */
         genericATdonorprob,          /*     generic prob. of AT donor */
         nongenericATdonorprob,       /* non-generic prob. of AT donor */
         genericAGacceptorprob,       /*     generic prob. of AG acceptor */
         nongenericAGacceptorprob,    /* non-generic prob. of AG acceptor */
         genericACacceptorprob,       /*     generic prob. of AC acceptor */
         nongenericACacceptorprob,    /* non-generic prob. of AC acceptor */
         genericothersplicesitep,     /*     generic prob. of other sp. sites */
         nongenericothersplicesitep;  /* non-generic prob. of other sp. sites */
  bool useU12intronmodel;             /* enable U12 intron model */
  GthFlt U12typedonorprob,            /* prob. of U12-type donor */
         U12typedonorprobonemismatch; /* prob. of U12-type donor with 1
                                         mismatch */

  /* the precomputed log values */
  GthFlt log_genericGTdonorprob,
         log1minus_genericGTdonorprob,
         log_nongenericGTdonorprob,
         log1minus_nongenericGTdonorprob,
         log_genericGCdonorprob,
         log1minus_genericGCdonorprob,
         log_nongenericGCdonorprob,
         log1minus_nongenericGCdonorprob,
         log_genericATdonorprob,
         log1minus_genericATdonorprob,
         log_nongenericATdonorprob,
         log1minus_nongenericATdonorprob,
         log_genericAGacceptorprob,
         log1minus_genericAGacceptorprob,
         log_nongenericAGacceptorprob,
         log1minus_nongenericAGacceptorprob,
         log_genericACacceptorprob,
         log1minus_genericACacceptorprob,
         log_nongenericACacceptorprob,
         log1minus_nongenericACacceptorprob,
         log_genericothersplicesitep,
         log1minus_genericothersplicesitep,
         log_nongenericothersplicesitep,
         log1minus_nongenericothersplicesitep,
         log_U12typedonorprob,
         log1minus_U12typedonorprob,
         log_U12typedonorprobonemismatch,
         log1minus_U12typedonorprobonemismatch;

  GthBSSMParam *bssm_param; /* contains bssm parameters or is NULL,
                              if generic parameters are used. */
};

#endif
