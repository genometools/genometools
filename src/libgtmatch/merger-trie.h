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

#ifndef MERGER_TRIE_H
#define MERGER_TRIE_H

#include "libgtcore/arraydef.h"
#include "seqpos-def.h"
#include "encseq-def.h"

/* XXX make this type opaque */

typedef struct
{
  unsigned int idx;
  Seqpos startpos;
#ifdef WITHTRIEIDENT
  uint64_t  ident;
#endif
} Suffixinfo;

typedef struct Mergertrienode
{
  Suffixinfo suffixinfo;
  struct Mergertrienode *firstchild,
                        *rightsibling,
                        *parent;
  Seqpos depth;
} Mergertrienode;

typedef Mergertrienode * Mergertrienodeptr;

DECLAREARRAYSTRUCT(Mergertrienodeptr);

typedef struct
{
  Encodedsequence *encseqptr;
  Readmode readmode;
} Encseqreadinfo;

typedef struct
{
  Encseqreadinfo *encseqreadinfo;
  Mergertrienode *nodetable,
           *root;
  Mergertrienodeptr *unusedMergertrienodes;
  unsigned int numofindexes,
               nextunused,
               allocatedMergertrienode,
               nextfreeMergertrienode;
} Mergertrierep;

void showmergertrie(const Mergertrierep *trierep,
                    const Uchar *characters);

void checkmergertrie(Mergertrierep *trierep,unsigned int numberofleaves,
                     unsigned int maxleafnum,Error *err);

void showallnoderelations(const Mergertrienode *node);

void insertsuffixintomergertrie(Mergertrierep *trierep,
                                Mergertrienode *node,
                                Suffixinfo *suffixinfo);

Mergertrienode *findsmallestnodeintrie(const Mergertrierep *trierep);

void deletesmallestpath(Mergertrienode *smallest,Mergertrierep *trierep);

void initmergertrienodetable(Mergertrierep *trierep,Seqpos numofsuffixes,
                       unsigned int numofindexes);

void freemergertrierep(Mergertrierep *trierep);

#endif
