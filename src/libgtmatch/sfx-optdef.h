#ifndef SFX_OPTDEF_H
#define SFX_OPTDEF_H

#include "libgtcore/str.h"

#define PREFIXLENGTH_AUTOMATIC 0

typedef struct
{
  const char **filenamelist;
  unsigned int numoffiles,
               numofparts,
               prefixlength;
  Str *str_indexname,
      *str_smap,
      *str_sat;
  bool isdna,
       isprotein,
       outsuftab,
       outbwttab,
       outtistab;
} Suffixeratoroptions;

#endif
