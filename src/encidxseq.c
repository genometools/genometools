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
#include <string.h>

#include <libgtmatch/sarr-def.h>
#include <libgtmatch/esa-map.pr>

#include "encidxseq.h"
#include "encidxseqpriv.h"

void
deleteEncIdxSeq(struct encIdxSeq *seq, Env *env)
{
  seq->classInfo->delete(seq, env);
}

#define verifyIntegrityErrRet(retval)                                   \
  do {                                                                  \
    switch(retval) {                                                    \
    case 1:                                                             \
      fprintf(stderr, "Comparision failed at position %llu"             \
              ", reference symbol: %u, symbol read: %u\n",              \
              (unsigned long long)pos, symOrig, symEnc);                \
      break;                                                            \
    case -1:                                                            \
      fprintf(stderr, "Read of symbol failed at position %llu\n",       \
              (unsigned long long)pos);                                 \
      break;                                                            \
    case 2:                                                             \
      fprintf(stderr, "At position %llu, rank operation yielded"        \
              " wrong count: "FormatSeqpos" expected "FormatSeqpos"\n", \
              (unsigned long long)pos, rankQueryResult, rankExpect);    \
      break;                                                            \
    }                                                                   \
    deleteEISHint(seqIdx, hint, env);                                   \
    freesuffixarray(&suffixArray, env);                                 \
    return retval;                                                      \
  } while(0)

/**
 * @param tickPrint if not zero, print a . every tickPrint symbols to
 * fp
 * @param fp descriptor to write progress marks to.
 * @return -1 on error, 0 on identity, >0 on inconsistency
 */
int
verifyIntegrity(struct encIdxSeq *seqIdx,
                Str *projectName, int tickPrint, FILE *fp, Env *env)
{
  Seqpos rankTable[UCHAR_MAX+1];
  FILE *bwtFP;
  off_t pos = 0;
  Suffixarray suffixArray;
  Symbol symOrig, symEnc;
  EISHint hint;
  Seqpos seqLastPos, rankQueryResult, rankExpect;
  /* two part process: enumerate all positions of original sequence
   * and verify that the query functions return correct values */
  if(streamsuffixarray(&suffixArray, &seqLastPos,
                       SARR_SUFTAB | SARR_BWTTAB, projectName, NULL, env))
    return -1;
  memset(rankTable, 0, sizeof(rankTable));
  bwtFP = suffixArray.bwttabstream.fp;
/*   pos = 1803218; */
/*   fseeko(bwtFP, pos, SEEK_SET); */
  hint = newEISHint(seqIdx, env);
  while((symOrig = getc(bwtFP)) != EOF)
  {
    /* TODO: complete once query functions are finished */
/*     fprintf(stderr, "pos: %llu\n", (unsigned long long)pos); */
    symEnc = EISGetSym(seqIdx, pos, hint, env);
    if(symEnc == EOF)
      verifyIntegrityErrRet(-1);
    if(symEnc != symOrig)
      verifyIntegrityErrRet(1);
    if((rankExpect = ++rankTable[symOrig])
       != (rankQueryResult = EISRank(seqIdx, symOrig, pos + 1, hint, env)))
      verifyIntegrityErrRet(2);
    ++pos;
    if(tickPrint && !(pos % tickPrint))
      putc('.', fp);
/*     if(pos > 2000) */
/*       break; */
  }
  if(tickPrint)
    putc('\n', fp);
  if(ferror(bwtFP))
    verifyIntegrityErrRet(-1);
  deleteEISHint(seqIdx, hint, env);
  freesuffixarray(&suffixArray, env);
  return 0;
}
