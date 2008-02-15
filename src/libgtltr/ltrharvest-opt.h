/*
  Copyright (c) 2007 David Ellinghaus <d.ellinghaus@ikmb.uni-kiel.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef LTRHARVEST_OPT_H
#define LTRHARVEST_OPT_H

#include <stdbool.h>
#include "libgtcore/error.h"
#include "libgtcore/str.h"
#include "libgtcore/strarray.h"
#include "libgtcore/symboldef.h"
#include "myxdrop.h"
#include "repeattypes.h"

typedef struct
{
  RepeatInfo repeatinfo;                  /* stores all repeats */
  ArrayLTRboundaries arrayLTRboundaries;  /* stores all predicted */
                                          /*   LTR elements */
  Seqpos *markpos;                        /* positions of SEPARATOR symbols */
                                          /* in encseq */

  Str *str_indexname;           /* name of the suffix array index */
  Str *str_fastaoutputfilename; /* name of the FASTA output file */
  Str *str_fastaoutputfilenameinnerregion;  /* name of the FASTA */
                                            /* file for the inner regions */
  Str *str_gff3filename;         /* name of the gff3 file */
  unsigned long minseedlength;   /* minimal exact seed */
  double similaritythreshold;    /* minimum similarity of LTRs */
  int xdropbelowscore;           /* xdropbelowscore */
  Arbitraryscores arbitscores;
  Motif motif;                   /* the start-, endmotiv of the LTRs, */
                                 /* by default: OFF */
  bool verbosemode;      /* show extra statements, by default: OFF */
  bool longoutput;       /* additionally shows motif and TSD infos */
                         /* by default: OFF */
  Str *str_overlaps;     /* string from argv */
  bool bestofoverlap;    /* take best prediction */
                         /* if overlap occurs, default */
  bool nooverlapallowed; /* overlapping predictions (not)allowed */
  bool fastaoutput;      /* by default no FASTA output */
  bool fastaoutputinnerregion;
  bool gff3output;       /* by default no gff3 output */

  unsigned int minlengthTSD,  /* minlength of TSD, default */
               maxlengthTSD;  /* maxlength of TSD, default */
                               /* by default no search for TSDs */
  Seqpos vicinityforcorrectboundaries; /* vicinity for search of TSD and
                                          motif */
} LTRharvestoptions;

void showuserdefinedoptionsandvalues(LTRharvestoptions *lo);

void printargsline(const char **argv, int argc);

int testmotifandencodemotif (Motif *motif, const Alphabet *alpha, Error *);

int ltrharvestoptions(LTRharvestoptions *lo, int argc,const char **argv,
                      Error *);

void wrapltrharvestoptions(LTRharvestoptions *lo);

#endif
