/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include "spacedef.h"
#include "sarr-def.h"
#include "trieins-def.h"
#include "encseq-def.h"
#include "alphadef.h"

#include "trieins.pr"
#include "esa-map.pr"

static void maketrie(Trierep *trierep,
                     /*@unused@*/ const Uchar *characters,
                     Seqpos len)
{
  Suffixinfo suffixinfo;

  suffixinfo.idx = 0;
#ifdef WITHTRIEIDENT
  suffixinfo.ident = 0;
#endif
  for (suffixinfo.startpos = 0;
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
  unsigned int numberofleaves = (unsigned int) seqlen+1;
  unsigned int maxleafnum = (unsigned int) seqlen;
#endif

  while (trierep->root != NULL && trierep->root->firstchild != NULL)
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
  bool haserr = false;
  Seqpos totallength;

  env_error_check(env);
  if (streamsuffixarray(&suffixarray,
                        &totallength,
                        SARR_ESQTAB,
                        indexname,
                        NULL,
                        env) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    Trierep trierep;
    const Uchar *characters;

    ALLOCASSIGNSPACE(trierep.encseqreadinfo,NULL,Encseqreadinfo,1);
    trierep.encseqreadinfo[0].encseqptr = suffixarray.encseq;
    trierep.encseqreadinfo[0].readmode = suffixarray.readmode;
    characters = getcharactersAlphabet(suffixarray.alpha);
    inittrienodetable(&trierep,totallength,(unsigned int) 1,env);
    maketrie(&trierep,characters,totallength);
    if (onlyins)
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
