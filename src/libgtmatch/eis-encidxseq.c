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

#include "libgtmatch/seqpos-def.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/esa-map.pr"

#include "libgtmatch/eis-encidxseq.h"
#include "libgtmatch/eis-encidxseqpriv.h"

void
deleteEncIdxSeq(EISeq *seq, Env *env)
{
  seq->classInfo->delete(seq, env);
}

const char *EISintegrityCheckResultStrings[] =
{
  "no error in bwt comparison.",
  "reading from the index produced an incorrect symbol.",
  "reading from the BWT reference file failed or delivered a symbol"
  " not in the alphabet (bwt file corrupt?) ",
  "loading/mapping of the suffix array project failed. \n"
  "(did you generate the BWT?)",
  "the rank operation delivered a wrong count"
};

#define verifyIntegrityErrRet(retval)                                   \
  do {                                                                  \
    switch (retval) {                                                   \
    case EIS_INTEGRITY_CHECK_INVALID_SYMBOL:                            \
      fprintf(stderr, "Comparision failed at position "FormatSeqpos     \
              ", reference symbol: %u, symbol read: %u\n",              \
              pos, symOrig, symEnc);                                    \
      break;                                                            \
    case EIS_INTEGRITY_CHECK_BWT_READ_ERROR:                            \
      fprintf(stderr, "Read of symbol failed at position "              \
              FormatSeqpos"\n", pos);                                   \
      break;                                                            \
    case EIS_INTEGRITY_CHECK_RANK_FAILED:                               \
      fprintf(stderr, "At position "FormatSeqpos                        \
              ", rank operation yielded  wrong count: "FormatSeqpos     \
              " expected "FormatSeqpos"\n",                             \
              pos, rankQueryResult, rankExpect);                        \
      break;                                                            \
    }                                                                   \
    EISPrintDiagsForPos(seqIdx, pos, stderr, hint, env);                \
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
extern enum integrityCheckResults
verifyIntegrity(EISeq *seqIdx, const Str *projectName, Seqpos skip,
                unsigned long tickPrint, FILE *fp, Env *env)
{
  Seqpos rankTable[UCHAR_MAX+1];
  FILE *bwtFP;
  Seqpos pos = 0;
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
    env_error_set(env, "Cannot load suffix array project with"
                  " demand for BWT file\n");
    freeverboseinfo(&verbosity, env);
    return EIS_INTEGRITY_CHECK_SA_LOAD_ERROR;
  }
  memset(rankTable, 0, sizeof (rankTable));
  bwtFP = suffixArray.bwttabstream.fp;
  hint = newEISHint(seqIdx, env);
  alphabet = EISGetAlphabet(seqIdx);
  if (skip > 0)
  {
    Seqpos len = EISLength(seqIdx);
    unsigned sym;
    if (skip >= len)
    {
      showverbose(verbosity, "Invalid skip request: %lld,"
                  " too large for sequence length: "FormatSeqpos,
                  (long long)skip, len);
      freeverboseinfo(&verbosity, env);
      return -1;
    }
    fseeko(bwtFP, skip, SEEK_SET);
    for (sym = 0; sym < UCHAR_MAX+1; ++sym)
      if (MRAEncSymbolHasValidMapping(alphabet, sym))
        rankTable[sym] = EISRank(seqIdx, sym, skip, hint, env);
    pos = skip;
  }
/*   EISPrintDiagsForPos(seqIdx, ((unsigned long)random())%EISLength(seqIdx), */
/*                       stderr, hint, env); */
  while ((symRead = getc(bwtFP)) != EOF)
  {
    symOrig = symRead;
    symEnc = EISGetSym(seqIdx, pos, hint, env);
    if (!MRAEncSymbolHasValidMapping(alphabet, symEnc))
      verifyIntegrityErrRet(EIS_INTEGRITY_CHECK_INVALID_SYMBOL);
    if (!MRAEncSymbolHasValidMapping(alphabet, symOrig))
      verifyIntegrityErrRet(EIS_INTEGRITY_CHECK_BWT_READ_ERROR);
    if (symEnc != symOrig)
      verifyIntegrityErrRet(EIS_INTEGRITY_CHECK_INVALID_SYMBOL);
    if ((rankExpect = ++rankTable[symOrig])
       != (rankQueryResult = EISRank(seqIdx, symOrig, pos + 1, hint, env)))
      verifyIntegrityErrRet(EIS_INTEGRITY_CHECK_RANK_FAILED);
    ++pos;
    if (tickPrint && !(pos % tickPrint))
      putc('.', fp);
  }
  if (tickPrint)
    putc('\n', fp);
  if (ferror(bwtFP))
    verifyIntegrityErrRet(EIS_INTEGRITY_CHECK_BWT_READ_ERROR);
  deleteEISHint(seqIdx, hint, env);
  freesuffixarray(&suffixArray, env);
  freeverboseinfo(&verbosity, env);
  return EIS_INTEGRITY_CHECK_NO_ERROR;
}
