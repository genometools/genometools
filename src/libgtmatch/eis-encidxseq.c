/*
  Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>

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
#include <string.h>

#include "libgtmatch/sarr-def.h"
#include "libgtmatch/esa-map.pr"

#include "libgtmatch/eis-encidxseq.h"
#include "libgtmatch/eis-encidxseqpriv.h"

void
deleteEncIdxSeq(EISeq *seq, Env *env)
{
  seq->classInfo->delete(seq, env);
}

#define verifyIntegrityErrRet(retval)                                   \
  do {                                                                  \
    switch (retval) {                                                    \
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
    freeverboseinfo(&verbosity, env);                                   \
    return retval;                                                      \
  } while (0)

/**
 * @param tickPrint if not zero, print a . every tickPrint symbols to
 * fp
 * @param fp descriptor to write progress marks to.
 * @return -1 on error, 0 on identity, >0 on inconsistency
 */
int
verifyIntegrity(EISeq *seqIdx, Str *projectName, int tickPrint,
                FILE *fp, Env *env)
{
  Seqpos rankTable[UCHAR_MAX+1];
  FILE *bwtFP;
  off_t pos = 0;
  Suffixarray suffixArray;
  Symbol symOrig, symEnc;
  EISHint hint;
  Seqpos seqLastPos, rankQueryResult, rankExpect;
  Verboseinfo *verbosity;
  const MRAEnc *alphabet;
  int symRead;
  verbosity = newverboseinfo(true, env);
  /* two part process: enumerate all positions of original sequence
   * and verify that the query functions return correct values */
  if (streamsuffixarray(&suffixArray, &seqLastPos,
                        SARR_BWTTAB, projectName, verbosity, env))
  {
    freeverboseinfo(&verbosity, env);
    return -1;
  }
  memset(rankTable, 0, sizeof (rankTable));
  bwtFP = suffixArray.bwttabstream.fp;
  alphabet = EISGetAlphabet(seqIdx);
/*   pos = 1803218; */
/*   fseeko(bwtFP, pos, SEEK_SET); */
  hint = newEISHint(seqIdx, env);
  while ((symRead = getc(bwtFP)) != EOF)
  {
    symOrig = symRead;
    /* TODO: complete once query functions are finished */
/*     fprintf(stderr, "pos: %llu\n", (unsigned long long)pos); */
    symEnc = EISGetSym(seqIdx, pos, hint, env);
    if (!MRAEncSymbolHasValidMapping(alphabet, symEnc))
      verifyIntegrityErrRet(-1);
    if (symEnc != symOrig)
      verifyIntegrityErrRet(1);
    if ((rankExpect = ++rankTable[symOrig])
       != (rankQueryResult = EISRank(seqIdx, symOrig, pos + 1, hint, env)))
      verifyIntegrityErrRet(2);
    ++pos;
    if (tickPrint && !(pos % tickPrint))
      putc('.', fp);
/*     if (pos > 2000) */
/*       break; */
  }
  if (tickPrint)
    putc('\n', fp);
  if (ferror(bwtFP))
    verifyIntegrityErrRet(-1);
  deleteEISHint(seqIdx, hint, env);
  freesuffixarray(&suffixArray, env);
  freeverboseinfo(&verbosity, env);
  return 0;
}
