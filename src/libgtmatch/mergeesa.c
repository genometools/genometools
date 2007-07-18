/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdio.h>
#include "sarr-def.h"
#include "emimergeesa.h"
#include "encseq-def.h"
#include "trieins-def.h"

#include "trieins.pr"
#include "sfxmap.pr"
#include "alphabet.pr"

 DECLAREREADFUNCTION(Seqpos);

 DECLAREREADFUNCTION(Uchar);

 DECLAREREADFUNCTION(Largelcpvalue);

static void fillandinsert(Trierep *trierep,
                          uint32_t idx,
                          Seqpos suftabvalue,
                          Trienode *node,
                          uint64_t ident)
{
  Suffixinfo sinfo;

  sinfo.idx = idx;
  sinfo.startpos = suftabvalue;
#ifdef DEBUG
  sinfo.ident = ident;
#endif
  insertsuffixintotrie(trierep,node,&sinfo);
}

static int inputthesequences(Seqpos *nextpostable,
                             Suffixarray *suffixarraytable,
                             const StrArray *indexnametab,
                             unsigned int demand,
                             Env *env)
{
  unsigned long idx;
  Str *indexname;
  
  for(idx=0; idx<strarray_size(indexnametab); idx++)
  {
    indexname = str_new_cstr(strarray_get(indexnametab,idx),env);
    if(streamsuffixarray(&suffixarraytable[idx],
                         demand,
                         indexname,
                         env) != 0)
    {
      str_delete(indexname,env);
      return -1;
    }
    str_delete(indexname,env);
    nextpostable[idx] = 0;
  }
  return 0;
}
