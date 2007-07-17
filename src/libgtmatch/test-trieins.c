/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "types.h"
#include "spacedef.h"
#include "sarr-def.h"
#include "trieins-def.h"
#include "encseq-def.h"

#include "trieins.pr"
#include "sfxmap.pr"
#include "alphabet.pr"

static void maketrie(Trierep *trierep,
                     /*@unused@*/ const Uchar *characters,
                     Seqpos len)
{
  Suffixinfo suffixinfo;

  suffixinfo.idx = 0;
#ifdef WITHIDENT
  suffixinfo.ident = 0;
#endif
  for(suffixinfo.startpos = 0; 
      suffixinfo.startpos <= len; 
      suffixinfo.startpos++)
  {
    insertsuffixintotrie(trierep,trierep->root,&suffixinfo);
#ifdef DEBUG
    DEBUG1(4,"\n\nstep%lu:\n",(Showuint) suffixinfo.startpos);
    showtrie(trierep,characters);
    suffixinfo.ident++;
#endif
  }
}

static void successivelydeletesmallest(Trierep *trierep,
                                       /*@unused@*/ Seqpos seqlen,
                                       /*@unused@*/ const Uchar *characters,
                                       /*@unused@*/ Env *env)
{
  Trienode *smallest;
#ifdef WITHTRIEIDENT
  Seqpos numberofleaves = seqlen;
#endif

  while(trierep->root != NULL && trierep->root->firstchild != NULL)
  {
    smallest = findsmallestnodeintrie(trierep);
    deletesmallestpath(smallest,trierep);
#ifdef WITHTRIEIDENT
    showtrie(trierep,characters);
    checktrie(trierep,numberofleaves,env);
    numberofleaves--;
#endif
  }
}

int test_trieins(bool onlyins,const Str *indexname,Env *env)
{
  Suffixarray suffixarray;
  Trierep trierep;
  bool haserr = true;
  Seqpos totallength;
  const Uchar *characters;

  if(streamsuffixarray(&suffixarray,
                       true,
                       false,
                       false,
                       false,
                       indexname,
                       env) != 0)
  {
    haserr = true;
  }
  if(!haserr)
  {
    trierep.encseqtable = &suffixarray.encseq;
    totallength = getencseqtotallength(suffixarray.encseq);
    characters = getcharactersAlphabet(suffixarray.alpha);
    inittrienodetable(&trierep,totallength,(uint32_t) 1,env);
    maketrie(&trierep,characters,totallength);
    if(onlyins)
    {
#ifdef WITHTRIEINS
      showtrie(&trierep,characters);
      checktrie(&trierep,totallength+1);
#endif
    } else
    {
      showallnoderelations(trierep.root);
      successivelydeletesmallest(&trierep,totallength,characters,env);
    }
    freetrierep(&trierep,env);
  }
  freesuffixarray(&suffixarray,env);
  return haserr ? -1 : 0;
}
