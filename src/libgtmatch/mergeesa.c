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
                          /*@unused@*/ uint64_t ident)
{
  Suffixinfo sinfo;

  sinfo.idx = idx;
  sinfo.startpos = suftabvalue;
#ifdef WITHTRIEIDENT
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
                               uint32_t numofindexes,
                               Env *env)
{
  uint32_t idx;
  Seqpos suftabvalue;
  int retval;

  for(idx=0; idx<numofindexes; idx++)
  {
    retval = readnextSeqposfromstream(&suftabvalue,
                                      &suffixarraytable[idx].suftabstream,
                                      env);
    if(retval < 0)
    {
      return -1; 
    }
    if(retval == 0)
    {
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
                                                       Seqpos lcpvalue,
                                                       Env *env)
{
  Trienode *tmp;

  for(tmp = smallest->parent; tmp != NULL; tmp = tmp->parent)
  {
    if(tmp->depth <= lcpvalue)
    {
      return tmp;
    }
  }
  env_error_set(env,"path does not contain a node of depth <= " FormatSeqpos,
                PRINTSeqposcast(lcpvalue));
  return NULL;
}

int stepdeleteandinsertothersuffixes(Emissionmergedesa *emmesa, Env *env)
{
  Trienode *tmpsmallestleaf, *tmplcpnode;
  Largelcpvalue tmpexception;
  Uchar tmpsmalllcpvalue;
  int retval;
  Seqpos tmpsuftabvalue,
         tmplcpvalue,
         tmplastbranchdepth;
  uint32_t tmpidx; 
  
  for(emmesa->buf.nextstoreidx = 0;
      emmesa->numofentries > 0 &&
      emmesa->buf.nextstoreidx < (uint32_t) SIZEOFMERGERESULTBUFFER;
      emmesa->buf.nextstoreidx++)
  {
    tmpsmallestleaf = findsmallestnodeintrie(&emmesa->trierep);
    tmplastbranchdepth = tmpsmallestleaf->parent->depth;
    tmpidx = tmpsmallestleaf->suffixinfo.idx;
    emmesa->buf.suftabstore[emmesa->buf.nextstoreidx].idx = tmpidx;
    emmesa->buf.suftabstore[emmesa->buf.nextstoreidx].startpos
      = tmpsmallestleaf->suffixinfo.startpos;
    if(emmesa->nextpostable[tmpidx] > 
       getencseqtotallength(emmesa->suffixarraytable[tmpidx].encseq))
    {
      deletesmallestpath(tmpsmallestleaf,&emmesa->trierep);
      emmesa->numofentries--;
    } else
    {
      retval = readnextUcharfromstream(&tmpsmalllcpvalue,
                                       &emmesa->suffixarraytable[tmpidx].
                                                lcptabstream,env);
      if(retval < 0)
      {
        return -1; 
      }
      if(retval == 0)
      {
        env_error_set(env,"file %s: line %d: unexpected end of file when "
                        "reading lcptab",__FILE__,__LINE__);
        return -2;
      }
      if(tmpsmalllcpvalue == (Uchar) UCHAR_MAX)
      {
        retval = readnextLargelcpvaluefromstream(
                               &tmpexception,
                               &emmesa->suffixarraytable[tmpidx].llvtabstream,
                               env);
        if(retval < 0)
        {
          return -3; 
        }
        if(retval == 0)
        {
          env_error_set(env,"file %s: line %d: unexpected end of file when "
                        "reading llvtab",__FILE__,__LINE__);
          return -4;
        }
        tmplcpvalue = tmpexception.value;
      } else
      {
        tmplcpvalue = (Seqpos) tmpsmalllcpvalue;
      }
      if(tmplcpvalue > tmplastbranchdepth)
      {
        tmplastbranchdepth = tmplcpvalue;
      }
      tmplcpnode = findlargestnodeleqlcpvalue(tmpsmallestleaf,tmplcpvalue,env);
      retval = readnextSeqposfromstream(&tmpsuftabvalue,
                                        &emmesa->suffixarraytable[tmpidx].
                                        suftabstream,
                                        env);
      if(retval < 0)
      {
        return -5;
      }
      if(retval == 0)
      {
        env_error_set(env,"file %s: line %d: unexpected end of file when "
                      "reading suftab",__FILE__,__LINE__);
        return -6;
      }
      emmesa->nextpostable[tmpidx]++;
      fillandinsert(&emmesa->trierep,
                    tmpidx,
                    tmpsuftabvalue,
                    tmplcpnode,
                    emmesa->ident++);
      tmpsmallestleaf = findsmallestnodeintrie(&emmesa->trierep);
      deletesmallestpath(tmpsmallestleaf,&emmesa->trierep);
    }
    if(emmesa->numofentries > 0)
    {
      emmesa->buf.lcptabstore[emmesa->buf.nextstoreidx]
        = tmplastbranchdepth;
      emmesa->buf.lastpage = false;
    } else
    {
      emmesa->buf.lastpage = true;
    }
  }
  return 0;
}

int initEmissionmergedesa(Emissionmergedesa *emmesa,
                          const StrArray *indexnametab,
                          unsigned int demand,
                          Env *env)
{
  uint32_t numofindexes;
  bool haserr = false;

  numofindexes = (uint32_t) strarray_size(indexnametab);
  emmesa->buf.nextaccessidx = emmesa->buf.nextstoreidx = 0;
  emmesa->numofindexes = numofindexes;
  emmesa->numofentries = numofindexes;
  emmesa->ident = (uint64_t) numofindexes;
  emmesa->trierep.encseqtable = NULL;
  ALLOCASSIGNSPACE(emmesa->suffixarraytable,NULL,Suffixarray,numofindexes);
  ALLOCASSIGNSPACE(emmesa->nextpostable,NULL,Seqpos,numofindexes);
  if(inputthesequences(emmesa->nextpostable,
                       emmesa->suffixarraytable,
                       indexnametab,
                       demand,
                       env) != 0)
  {
    haserr = true;
    return -1;
  }
  if(!haserr && numofindexes > (uint32_t) 1)
  {
    uint32_t idx;

    ALLOCASSIGNSPACE(emmesa->trierep.encseqtable,NULL,Encodedsequence *,
                     numofindexes);
    for(idx = 0; idx < numofindexes; idx++)
    {
      emmesa->trierep.encseqtable[idx] = emmesa->suffixarraytable[idx].encseq;
    }
    inittrienodetable(&emmesa->trierep,numofindexes,numofindexes,env);
    if(insertfirstsuffixes(&emmesa->trierep,
                           emmesa->nextpostable,
                           emmesa->suffixarraytable,
                           numofindexes,
                           env) != 0)
    {
      FREESPACE(emmesa->trierep.encseqtable);
      haserr = true;
    }
  }
  if(haserr)
  {
    FREESPACE(emmesa->suffixarraytable);
    FREESPACE(emmesa->nextpostable);
  }
  return haserr ? -1 : 0;
}

void wraptEmissionmergedesa(Emissionmergedesa *emmesa,Env *env)
{
  uint32_t idx;

  for(idx = 0; idx < emmesa->numofindexes; idx++)
  {
    freesuffixarray(emmesa->suffixarraytable + idx,env);
  }
  FREESPACE(emmesa->suffixarraytable);
  FREESPACE(emmesa->trierep.encseqtable);
  if(emmesa->numofindexes > (uint32_t) 1)
  {
    freetrierep(&emmesa->trierep,env);
  }
  FREESPACE(emmesa->nextpostable);
}
