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

static int insertfirstsuffixes(Trierep *trierep,
                               Seqpos *nextpostable,
                               Suffixarray *suffixarraytable,
                               uint32_t numofindexes)
{
  uint32_t idx;
  Seqpos suftabvalue;
  int retval;

  for(idx=0; idx<numofindexes; idx++)
  {
    retval = readnextSeqpostfromstream(&suftabvalue,
                                       &suffixarraytab[idx].suftabstream,
                                       env);
    if(retval < 0)
    {
      return -1; 
    }
    if(retval == 0)
    {
      ERROREOF("suftab");
      env_error_set(env,"file %s: line %d: unexpected end of file when "
                        "reading suftab",__FILE__,__LINE__);
      return -2;
    }
    nextpostable[idx]++;
    fillandinsert(trierep,
                  idx,
                  suftabvalue,
                  trierep->root,
                  (uint64_t) idx);
  }
  return 0;
}

/*@null@*/ static Trienode *findlargestnodeleqlcpvalue(Trienode *smallest,
                                                       Seqpos lcpvalue)
{
  Trienode *tmp;

  for(tmp = smallest->parent; tmp != NULL; tmp = tmp->parent)
  {
    if(tmp->depth <= lcpvalue)
    {
      return tmp;
    }
  }
  ERROR1("path does not contain a node of depth <= %lu",
          (Showuint) lcpvalue);
  return NULL;
}
