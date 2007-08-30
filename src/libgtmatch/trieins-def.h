/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef TRIEINS_DEF_H
#define TRIEINS_DEF_H

#include "seqpos-def.h"
#include "arraydef.h"
#include "encseq-def.h"

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

#endif
