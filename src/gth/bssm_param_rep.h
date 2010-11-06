/*
  Copyright (c) 2003-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2005 Michael E Sparks <mespar1@iastate.edu>
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

#ifndef BSSM_PARAM_REP_H
#define BSSM_PARAM_REP_H

#include "gth/bssm_param.h"

#define WINSIZE         100 /* 50nt on the left and right of GT|AG */
#define HYPOTHESIS7     7   /* T1, T2, T0, F1, F2, F0 and Fi */
#define HYPOTHESIS2     2   /* true (Tr) or false (fs) sites */

/* The following structures define the tables for the BSSM parameters. */
typedef GthFlt Hypo2table[HYPOTHESIS2][WINSIZE+2][4][4];
typedef GthFlt Hypo7table[HYPOTHESIS7][WINSIZE+2][4][4];

/* The version of the Bssmparam structure */
#define BSSMPARAMVERSION        2

typedef struct {
  unsigned long hypothesis_num;  /* number of hypothesis, either HYPOTHESIS2 or
                                    HYPOTHESIS7 */
  unsigned long window_size_left,
                window_size_right;
  union {
    Hypo2table hypo2table;
    Hypo7table hypo7table;
  } hypotables;
} GthBSSMModel;

struct GthBSSMParam{
  unsigned char version_num;   /* contains version number of the BSSM parameter
                                  structure */
  bool gt_donor_model_set,     /* use GT donor site model */
       gc_donor_model_set,     /* use GC donor site model */
       ag_acceptor_model_set;  /* use AG acceptor site model */
  GthBSSMModel gt_donor_model,
               gc_donor_model,
               ag_acceptor_model;
};

#endif
