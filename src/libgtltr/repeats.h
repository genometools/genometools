/*
  Copyright (C) 2007 David Ellinghaus <dellinghaus@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef REPEATS_H
#define REPEATS_H 

#include "libgtmatch/sarr-def.h"
#include "repeattypes.h"
#include "ltrharvest-opt.h"

void showonstdout (char *s);

void showrepeats (RepeatInfo * repeatinfo,
                  unsigned long seedminlength);

int simpleexactselfmatchoutput(/*@unused@*/ void *info,
    Seqpos len, Seqpos pos1, Seqpos pos2);

int simpleexactselfmatchstore(LTRharvestoptions *lo, Seqpos len, 
                              Seqpos pos1, Seqpos pos2); 

int subsimpleexactselfmatchstore (SubRepeatInfo * info, Seqpos len,
                                  Seqpos pos1, Seqpos pos2);
#endif
