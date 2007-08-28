/*
  Copyright (C) 2007 by David Ellinghaus <dellinghaus@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef LTRHARVEST_OPTDEF
#define LTRHARVEST_OPTDEF

#include <stdbool.h>

/* Prof. Kurtz */
//#include "virtualdef.h" //muss ersetzt werden

/* Gremme */
#include "libgtcore/str.h"
#include "libgtcore/strarray.h"

/* Ellinghaus */
#include "myxdrop.h"
#include "repeattypes.h"
#include "searchforLTRs.h"

typedef struct
{
  //Virtualtree virtualtree; // ersetzen
  //unsigned int demand;
  RepeatInfo repeatinfo;                  /* stores all repeats */
  ArrayLTRboundaries arrayLTRboundaries;  /* stores all predicted 
                                             LTR elements */
  Str *str_indexname;           /* name of the suffix array index */
  Str *str_fastaoutputfilename; /* name of the FASTA output file */
  Str *str_fastaoutputfilenameinnerregion;  /* name of the FASTA
                                               file for the inner regions */
  Str *str_gff3filename;        /* name of the gff3 file*/
  unsigned long minseedlength;   /* minimal exact seed */
  double similaritythreshold;    /* minimum similarity of LTRs */
  unsigned int xdropbelowscore;  /* xdropbelowscore */
  Arbitraryscores arbitscores;
  Motif motif;                   /* the start-, endmotiv of the LTRs,
                                    by default: OFF*/
  bool verbosemode;      /* show extra statements, by default: OFF */
  bool longoutput;       /* additionally shows motif and TSD infos */
			 /* by default: OFF */
  Str *str_overlaps;     /* string from argv */
  bool bestofoverlap;    /* take best prediction 
                                    if overlap occurs, default */
  bool nooverlapallowed; /* overlapping predictions (not)allowed*/
  bool fastaoutput;      /* by default no FASTA output */        
  bool fastaoutputinnerregion;
  bool gff3output;       /* by default no gff3 output */
  
  unsigned long minlengthTSD,  /* minlength of TSD, default */
                maxlengthTSD;  /* maxlength of TSD, default */
                               /* by default no search for TSDs */
  unsigned int vicinityforcorrectboundaries; /* vicinity for search of TSD
					        and motif */
} LTRharvestoptions;

int ltrharvestoptions(LTRharvestoptions *lo, int argc, 
                        const char **argv, Env *env);
void wrapltrharvestoptions(LTRharvestoptions *lo,Env *env);
#endif
