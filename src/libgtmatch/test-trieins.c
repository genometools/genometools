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
#ifdef WITHTRIEIDENT
  suffixinfo.ident = 0;
#endif
  for(suffixinfo.startpos = 0; 
      suffixinfo.startpos <= len; 
      suffixinfo.startpos++)
  {
    insertsuffixintotrie(trierep,trierep->root,&suffixinfo);
#ifdef WITHTRIEIDENT
#ifdef WITHTRIESHOW
    showtrie(trierep,characters);
#endif
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
  uint32_t numberofleaves = (uint32_t) seqlen+1;
  uint32_t maxleafnum = (uint32_t) seqlen;
#endif

  while(trierep->root != NULL && trierep->root->firstchild != NULL)
  {
    smallest = findsmallestnodeintrie(trierep);
    deletesmallestpath(smallest,trierep);
#ifdef WITHTRIEIDENT
#ifdef WITHTRIESHOW
    showtrie(trierep,characters);
#endif
    numberofleaves--;
    checktrie(trierep,numberofleaves,maxleafnum,env);
#endif
  }
}

int test_trieins(bool onlyins,const Str *indexname,Env *env)
{
  Suffixarray suffixarray;
  Trierep trierep;
  bool haserr = false;
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
#ifdef WITHTRIEIDENT
#ifdef WITHTRIESHOW
      showtrie(&trierep,characters);
#endif
      checktrie(&trierep,totallength+1,totallength,env);
#endif
    } else
    {
#ifdef WITHTRIEIDENT
#ifdef WITHTRIESHOW
      showallnoderelations(trierep.root);
#endif
#endif
      successivelydeletesmallest(&trierep,totallength,characters,env);
    }
    freetrierep(&trierep,env);
  }
  freesuffixarray(&suffixarray,env);
  return haserr ? -1 : 0;
}
