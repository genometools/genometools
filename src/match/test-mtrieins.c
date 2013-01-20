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

#include "core/alphabet.h"
#include "core/unused_api.h"
#include "core/encseq.h"
#include "core/ma_api.h"
#include "sarr-def.h"
#include "merger-trie.h"
#include "esa-map.h"

static void maketrie(Mergertrierep *trierep,
                     GT_UNUSED const GtUchar *characters,
                     unsigned long len)
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
    gt_mergertrie_insertsuffix(trierep,trierep->root,&suffixinfo);
#ifdef WITHTRIEIDENT
#ifdef WITHTRIESHOW
    showtrie(trierep,characters);
#endif
    suffixinfo.ident++;
#endif
  }
}

static void successivelydeletesmallest(Mergertrierep *trierep,
                                       GT_UNUSED unsigned long seqlen,
                                       GT_UNUSED const GtUchar *characters,
                                       GT_UNUSED GtError *err)
{
  Mergertrienode *smallest;
#ifdef WITHTRIEIDENT
  unsigned int numberofleaves = (unsigned int) seqlen+1;
  unsigned int maxleafnum = (unsigned int) seqlen;
#endif

  while (trierep->root != NULL && trierep->root->firstchild != NULL)
  {
    smallest = gt_mergertrie_findsmallestnode(trierep);
    gt_mergertrie_deletesmallestpath(smallest,trierep);
#ifdef WITHTRIEIDENT
#ifdef WITHTRIESHOW
    showtrie(trierep,characters);
#endif
    numberofleaves--;
    checktrie(trierep,numberofleaves,maxleafnum,err);
#endif
  }
}

int gt_test_trieins(bool onlyins,const char *indexname,GtError *err)
{
  Suffixarray suffixarray;
  bool haserr = false;
  unsigned long totallength = 0;

  gt_error_check(err);
  if (streamsuffixarray(&suffixarray,
                        SARR_ESQTAB,
                        indexname,
                        NULL,
                        err) != 0)
  {
    haserr = true;
  } else
  {
    totallength = gt_encseq_total_length(suffixarray.encseq);
  }
  if (!haserr)
  {
    Mergertrierep trierep;
    const GtUchar *characters;

    trierep.encseqreadinfo = gt_malloc(sizeof *trierep.encseqreadinfo);
    trierep.encseqreadinfo->encseqptr = suffixarray.encseq;
    trierep.encseqreadinfo->readmode = suffixarray.readmode;
    characters
      = gt_alphabet_characters(gt_encseq_alphabet(suffixarray.encseq));
    gt_mergertrie_initnodetable(&trierep,totallength,1U);
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
    gt_mergertrie_delete(&trierep);
  }
  gt_freesuffixarray(&suffixarray);
  return haserr ? -1 : 0;
}
