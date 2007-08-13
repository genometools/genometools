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

#include "encidxseq.h"
#include "encidxseqpriv.h"

staticifinline inline Seqpos
EISLength(struct encIdxSeq *seq)
{
  return seq->seqLen;
}

staticifinline inline Symbol
EISGetSym(struct encIdxSeq *seq, Seqpos pos, union EISHint *hint, Env *env)
{
  return seq->classInfo->get(seq, pos, hint, env);
}

staticifinline inline Seqpos
EISRank(struct encIdxSeq *seq, Symbol sym, Seqpos pos, union EISHint *hint,
        Env *env)
{
  return seq->classInfo->rank(seq, sym, pos, hint, env);
}

staticifinline inline EISHint
newEISHint(struct encIdxSeq *seq, Env *env)
{
  return seq->classInfo->newHint(seq, env);
}

staticifinline inline void
deleteEISHint(struct encIdxSeq *seq, EISHint hint, Env *env)
{
  return seq->classInfo->deleteHint(seq, hint, env);
}


