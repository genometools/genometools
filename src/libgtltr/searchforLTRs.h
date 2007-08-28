/* Copyright (C) 2007 by David Ellinghaus <dellinghaus@zbh.uni-hamburg.de> */

#ifndef SEARCHFORLTRS_H
#define SEARCHFORLTRS_H

//#include "types.h"
#include <libgtcore/str.h> 

/*
 The datatype Motif stores information about the specified motif.
 */
typedef struct
{
  Str *str_motif;
  unsigned char firstleft,     /* first character of left motif instance */
        secondleft,    /* second character of left motif instance */ 
	firstright,    /* first character of right motif instance */ 
	secondright;   /* second character of right motif instance */ 
  unsigned int allowedmismatches; /* number of allowed mismatches in the four 
                             character motif*/
} Motif;
#endif
