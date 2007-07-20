/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef SFX_OPTDEF_H
#define SFX_OPTDEF_H

#include "libgtcore/str.h"
#include "libgtcore/strarray.h"

#define PREFIXLENGTH_AUTOMATIC 0

typedef struct
{
  unsigned int numofparts,
               prefixlength;
  Str *str_indexname,
      *str_smap,
      *str_sat;
  StrArray *filenametab;
  bool isdna,
       isprotein,
       isplain,
       outtistab,
       outsuftab,
       outlcptab,
       outbwttab;
} Suffixeratoroptions;

#endif
