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

#include "libgtcore/unused.h"
#include "spacedef.h"
#include "sarr-def.h"
#include "merger-trie.h"
#include "encseq-def.h"
#include "alphadef.h"

#include "esa-map.pr"

static void maketrie(Mergertrierep *trierep,
                     UNUSED const Uchar *characters,
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
    insertsuffixintomergertrie(trierep,trierep->root,&suffixinfo);
#ifdef WITHTRIEIDENT
#ifdef WITHTRIESHOW
    showtrie(trierep,characters);
#endif
    suffixinfo.ident++;
#endif
  }
}

static void successivelydeletesmallest(Mergertrierep *trierep,
                                       UNUSED Seqpos seqlen,
                                       UNUSED const Uchar *characters,
                                       UNUSED Error *err)
{
  Mergertrienode *smallest;
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
    checktrie(trierep,numberofleaves,maxleafnum,err);
#endif
  }
}

int test_trieins(bool onlyins,const Str *indexname,Error *err)
{
  Suffixarray suffixarray;
  bool haserr = false;
  Seqpos totallength;

  error_check(err);
  if (streamsuffixarray(&suffixarray,
                        &totallength,
                        SARR_ESQTAB,
                        indexname,
                        NULL,
                        err) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    Mergertrierep trierep;
    const Uchar *characters;

    ALLOCASSIGNSPACE(trierep.encseqreadinfo,NULL,Encseqreadinfo,1);
    trierep.encseqreadinfo[0].encseqptr = suffixarray.encseq;
    trierep.encseqreadinfo[0].readmode = suffixarray.readmode;
    characters = getcharactersAlphabet(suffixarray.alpha);
    initmergertrienodetable(&trierep,totallength,1U);
    maketrie(&trierep,characters,totallength);
    if (onlyins)
    {
#ifdef WITHTRIEIDENT
#ifdef WITHTRIESHOW
      showtrie(&trierep,characters);
#endif
      checktrie(&trierep,totallength+1,totallength,err);
#endif
    } else
    {
#ifdef WITHTRIEIDENT
#ifdef WITHTRIESHOW
      showallnoderelations(trierep.root);
#endif
#endif
      successivelydeletesmallest(&trierep,totallength,characters,err);
    }
    freemergertrierep(&trierep);
  }
  freesuffixarray(&suffixarray);
  return haserr ? -1 : 0;
}
