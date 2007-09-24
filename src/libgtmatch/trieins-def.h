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

#ifndef TRIEINS_DEF_H
#define TRIEINS_DEF_H

#include "seqpos-def.h"
#include "arraydef.h"
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

typedef struct _Trienode
{
  Suffixinfo suffixinfo;
  struct _Trienode *firstchild,
                   *rightsibling,
                   *parent;
  Seqpos depth;
} Trienode;

typedef Trienode * Trienodeptr;

DECLAREARRAYSTRUCT(Trienodeptr);

typedef struct
{
  Encodedsequence *encseqptr;
  Readmode readmode;
} Encseqreadinfo;

typedef struct
{
  Encseqreadinfo *encseqreadinfo;
  Trienode *nodetable,
           *root;
  Trienodeptr *unusedTrienodes;
  unsigned int numofindexes,
               nextunused,
               allocatedTrienode,
               nextfreeTrienode;
} Trierep;

void showtrie(const Trierep *trierep,
              const Uchar *characters);

void checktrie(Trierep *trierep,unsigned int numberofleaves,
               unsigned int maxleafnum,Env *env);

void showallnoderelations(const Trienode *node);

void insertsuffixintotrie(Trierep *trierep,
                          Trienode *node,
                          Suffixinfo *suffixinfo);

Trienode *findsmallestnodeintrie(const Trierep *trierep);

void deletesmallestpath(Trienode *smallest,Trierep *trierep);

void inittrienodetable(Trierep *trierep,Seqpos numofsuffixes,
                       unsigned int numofindexes,Env *env);

void freetrierep(Trierep *trierep,Env *env);

#endif
