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
#ifndef ENCIDXSEQ_H_INCLUDED
#define ENCIDXSEQ_H_INCLUDED

/**
 * \file encidxseq.h
 * Interface definitions for encoded indexed sequences.
 */

#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif /* HAVE_STDINT_H */
#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#endif /* HAVE_INTTYPES_H */

#include <libgtcore/env.h>
#include <libgtcore/str.h>
#include <libgtmatch/seqpos-def.h>

#include "mrangealphabet.h"

/* typedef int_fast64_t Seqpos; */
enum rangeStoreMode {
  DIRECT_SYM_ENCODE,
  BLOCK_COMPOSITION_INCLUDE,
  REGIONS_LIST,
};

typedef union EISHint *EISHint;

extern struct encIdxSeq *
newBlockEncIdxSeq(enum rangeStoreMode modes[],
                  Str *projectName,
                  unsigned blockSize, Env *env);
struct encIdxSeq *
loadBlockEncIdxSeq(Str *projectName, Env *env);

extern int
searchBlock2IndexPair(const struct encIdxSeq *seqIdx,
                      const Symbol *block,
                      size_t idxOutput[2], Env *env);

extern void
deleteEncIdxSeq(struct encIdxSeq *seq, Env *env);

staticifinline inline Seqpos
EISRank(struct encIdxSeq *seq, Symbol sym, Seqpos pos, union EISHint *hint,
        Env *env);

extern Seqpos
EISSelect(struct encIdxSeq *seq, Symbol sym, Seqpos count);

staticifinline inline Seqpos
EISLength(struct encIdxSeq *seq);

staticifinline inline Symbol
EISGetSym(struct encIdxSeq *seq, Seqpos pos, EISHint hint, Env *env);

staticifinline inline EISHint
newEISHint(struct encIdxSeq *seq, Env *env);

staticifinline inline void
deleteEISHint(struct encIdxSeq *seq, EISHint hint, Env *env);

extern int
verifyIntegrity(struct encIdxSeq *seqIdx,
                Str *projectName, int tickPrint, FILE *fp, Env *env);

#ifdef HAVE_WORKING_INLINE
#include "../encidxseqsimpleop.c"
#endif

#endif /* ENCIDXSEQ_H_INCLUDED */
