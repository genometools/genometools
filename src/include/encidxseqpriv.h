/*
** Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>
**  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
**  
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**  
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
**  
*/
#ifndef ENCIDXSEQPRIV_H_INCLUDED
#define ENCIDXSEQPRIV_H_INCLUDED

#include "seqranges.h"
#include "encidxseq.h"

struct encIdxSeqClass
{
  void (*delete)(struct encIdxSeq *seq, Env *env);
  Seqpos (*rank)(struct encIdxSeq *seq, Symbol sym, Seqpos pos,
                 union EISHint *hint, Env *env);
  Seqpos (*select)(struct encIdxSeq *seq, Symbol sym, Seqpos count,
                   union EISHint *hint, Env *env);
  Symbol (*get)(struct encIdxSeq *seq, Seqpos pos, EISHint hint,
                Env *env);
  union EISHint *(*newHint)(struct encIdxSeq *seq, Env *env);
  void (*deleteHint)(struct encIdxSeq *seq, EISHint hint, Env *env);

};

struct encIdxSeq
{
  const struct encIdxSeqClass *classInfo;
  Seqpos seqLen;
};

struct seqCache
{
  size_t numEntries;
  Seqpos *cachedPos;
  void **entriesPtr;
  void *entries;
};


struct blockEncIdxSeqHint
{
  struct seqCache sBlockCache;
  seqRangeListSearchHint rangeHint;
};

union EISHint
{
  struct blockEncIdxSeqHint bcHint;
};

#endif /* ENCIDXSEQPRIV_H_INCLUDED */
