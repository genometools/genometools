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

#ifndef BSSM_PARAM_H
#define BSSM_PARAM_H

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "core/error.h"
#include "core/file.h"

typedef struct GthBSSMParam GthBSSMParam;

typedef float  GthFlt; /* probability type for low precision */
typedef double GthDbl; /* probability type for high precision */

typedef enum {
  GT_DONOR_TYPE = 0,
  GC_DONOR_TYPE,
  AG_ACCEPTOR_TYPE,
  NUMOFTYPES
} Termtype;

#define GTH_UNDEF_GTHFLT  FLT_MAX
#define GTH_UNDEF_GTHDBL  DBL_MAX

#define BSSMFILEENDING          "bssm"

GthBSSMParam* gth_bssm_param_new(void);
/* Read a bssm parameter file on the path $BSSMDIR. */
GthBSSMParam* gth_bssm_param_load(const char *filename, GtError*);
GthBSSMParam* gth_bssm_param_extract(unsigned long speciesnum, GtError*);
void          gth_bssm_param_delete(GthBSSMParam*);
/* Save the data contained in <bssm_param> to a file named <filename>. */
int           gth_bssm_param_save(GthBSSMParam*, const char *filename,
                                  GtError*);
/* Returns <true>, if <bssm_param> is a seven-class model. <false> otherwise. */
bool          gth_bssm_param_is_seven_class(const GthBSSMParam *bssm_param);
/* Prints the contents of the bssm parameterization <bssm_param> to <outfp>. */
void          gth_bssm_param_echo(const GthBSSMParam *bssm_param, FILE *outfp);
void          gth_bssm_param_show_info(const GthBSSMParam*, GtFile *outfp);
int           gth_bssm_param_parameterize(GthBSSMParam*, const char *path,
                                          Termtype, bool gzip, GtError*);

#endif
