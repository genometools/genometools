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
#ifndef BWTSEQ_H_INCLUDED
#define BWTSEQ_H_INCLUDED

#include <libgtcore/env.h>
#include <libgtcore/str.h>
#include <libgtmatch/seqpos-def.h>

#include "mrangealphabet.h"

/* FIXME:
 * - implement other index types
 */

enum seqBaseEncoding {
  BWT_BASETYPE_AUTOSELECT,
  BWT_ON_RLE,
  BWT_ON_BLOCK_ENC,
  BWT_ON_WAVELET_TREE_ENC,
};

union bwtSeqParam
{
  struct 
  {
    unsigned blockSize;
  } blockEncParams;
/*   struct  */
/*   { */
/*   } RLEParams; */
/*   struct  */
/*   { */
/*   } waveletTreeParams; */
};


typedef struct BWTSeq BWTSeq;

extern BWTSeq *
newBWTSeq(enum seqBaseEncoding baseType, union bwtSeqParam *extraParams,
          const Str *projectName, Env *env);

extern void
deleteBWTSeq(BWTSeq *bwtseq, Env *env);

extern int
BWTSeqHasLocateInformation(const BWTSeq *bwtseq);

staticifinline inline const MRAEnc *
BWTSeqGetAlphabet(const BWTSeq *seq);

staticifinline inline Seqpos
BWTSeqLength(const BWTSeq *seq);

extern Seqpos
BWTSeqMatchCount(const BWTSeq *bwtseq, Symbol *query, size_t queryLen,
                 Env *env);

extern struct BWTSeqExactMatchesIterator *
newEMIterator(const BWTSeq *bwtSeq, Symbol *query, size_t queryLen, Env *env);

extern void
deleteEMIterator(struct BWTSeqExactMatchesIterator *iter, Env *env);

struct MatchData
{
  const char *dbFile;
  Seqpos sfxArrayValue, dbFilePos;
};

extern struct MatchData *
EMIGetNextMatch(struct BWTSeqExactMatchesIterator *iter, const BWTSeq *bwtSeq,
                Env *env);

extern Seqpos
EMINumMatchesTotal(const struct BWTSeqExactMatchesIterator *iter);

extern Seqpos
EMINumMatchesLeft(const struct BWTSeqExactMatchesIterator *iter);

#ifdef HAVE_WORKING_INLINE
#include "../bwtseqsimpleop.c"
#endif

#endif /* BWTSEQ_H_INCLUDED */
