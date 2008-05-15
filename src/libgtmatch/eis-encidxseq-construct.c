/*
  Copyright (C) 2008 Thomas Jahns <Thomas.Jahns@gmx.net>

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

#include "libgtmatch/eis-blockcomp-construct.h"
#include "libgtmatch/eis-encidxseq-construct.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/esa-map.pr"

static EISeq *
createEncIdxSeqFromSASeqSrc(SASeqSrc *src,
                            const Str *projectName,
                            const struct seqBaseParam *params,
                            size_t numExtHeaders, const uint16_t *headerIDs,
                            const uint32_t *extHeaderSizes,
                            headerWriteFunc *extHeaderCallbacks,
                            void **headerCBData, bitInsertFunc biFunc,
                            BitOffset cwExtBitsPerPos,
                            varExtBitsEstimator biVarBits, void *cbState,
                            Error *err);

extern EISeq *
createEncIdxSeq(const Str *projectName,
                const struct seqBaseParam *params,
                size_t numExtHeaders, const uint16_t *headerIDs,
                const uint32_t *extHeaderSizes,
                headerWriteFunc *extHeaderCallbacks, void **headerCBData,
                bitInsertFunc biFunc, BitOffset cwExtBitsPerPos,
                varExtBitsEstimator biVarBits, void *cbState,
                Verboseinfo *verbosity, Error *err)
{
  Suffixarray suffixArray;
  struct encIdxSeq *newSeqIdx;
  Seqpos length;
  assert(projectName);
  /* map and interpret index project file */
  if (streamsuffixarray(&suffixArray, &length,
                       SARR_SUFTAB | SARR_BWTTAB, projectName, verbosity, err))
    return NULL;
  ++length;
  newSeqIdx = createEncIdxSeqFromSA(&suffixArray, length,
                                      projectName, params,
                                      numExtHeaders, headerIDs,
                                      extHeaderSizes, extHeaderCallbacks,
                                      headerCBData, biFunc, cwExtBitsPerPos,
                                      biVarBits, cbState, err);
  freesuffixarray(&suffixArray);
  return newSeqIdx;
}

extern EISeq *
createEncIdxSeqFromSA(Suffixarray *sa, Seqpos totalLen,
                      const Str *projectName,
                      const struct seqBaseParam *params,
                      size_t numExtHeaders, const uint16_t *headerIDs,
                      const uint32_t *extHeaderSizes,
                      headerWriteFunc *extHeaderCallbacks,
                      void **headerCBData,
                      bitInsertFunc biFunc, BitOffset cwExtBitsPerPos,
                      varExtBitsEstimator biVarBits, void *cbState,
                      Error *err)
{
  struct encIdxSeq *newSeqIdx;
  SuffixarrayFileInterface sai;
  assert(sa && projectName && err);
  initSuffixarrayFileInterface(&sai, totalLen, sa);
  newSeqIdx = createEncIdxSeqFromSAI(
    &sai, projectName, params, numExtHeaders, headerIDs,
    extHeaderSizes, extHeaderCallbacks, headerCBData, biFunc, cwExtBitsPerPos,
    biVarBits, cbState, err);
  destructSuffixarrayFileInterface(&sai);
  return newSeqIdx;
}

extern EISeq *
createEncIdxSeqFromSAI(SuffixarrayFileInterface *sai,
                       const Str *projectName,
                       const struct seqBaseParam *params,
                       size_t numExtHeaders, const uint16_t *headerIDs,
                       const uint32_t *extHeaderSizes,
                       headerWriteFunc *extHeaderCallbacks,
                       void **headerCBData,
                       bitInsertFunc biFunc, BitOffset cwExtBitsPerPos,
                       varExtBitsEstimator biVarBits, void *cbState,
                       Error *err)
{
  assert(sai && projectName && err);
  return createEncIdxSeqFromSASeqSrc(
    SAI2SASS(sai), projectName, params,
    numExtHeaders, headerIDs, extHeaderSizes,
    extHeaderCallbacks, headerCBData, biFunc,
    cwExtBitsPerPos, biVarBits, cbState, err);
}

extern EISeq *
createEncIdxSeqFromSfxI(sfxInterface *sfxi,
                        const Str *projectName,
                        const struct seqBaseParam *params,
                        size_t numExtHeaders, const uint16_t *headerIDs,
                        const uint32_t *extHeaderSizes,
                        headerWriteFunc *extHeaderCallbacks,
                        void **headerCBData,
                        bitInsertFunc biFunc, BitOffset cwExtBitsPerPos,
                        varExtBitsEstimator biVarBits, void *cbState,
                        Error *err)
{
  assert(sfxi && projectName && err);
  return createEncIdxSeqFromSASeqSrc(
    SfxI2SASS(sfxi), projectName, params,
    numExtHeaders, headerIDs, extHeaderSizes,
    extHeaderCallbacks, headerCBData, biFunc,
    cwExtBitsPerPos, biVarBits, cbState, err);
}

static EISeq *
createEncIdxSeqFromSASeqSrc(SASeqSrc *src,
                            const Str *projectName,
                            const struct seqBaseParam *params,
                            size_t numExtHeaders, const uint16_t *headerIDs,
                            const uint32_t *extHeaderSizes,
                            headerWriteFunc *extHeaderCallbacks,
                            void **headerCBData, bitInsertFunc biFunc,
                            BitOffset cwExtBitsPerPos,
                            varExtBitsEstimator biVarBits, void *cbState,
                            Error *err)
{
  SeqDataReader readSfxBWTSym;
  MRAEnc *alphabet;
  struct encIdxSeq *newSeqIdx;
  if (!SDRIsValid(readSfxBWTSym
                  = SASSCreateReader(src, SFX_REQUEST_BWTTAB)))
    return NULL;
  alphabet = SASSNewMRAEnc(src);
  newSeqIdx = createEncIdxSeqGen(SASSGetLength(src), projectName,
                                 alphabet, SASSGetSeqStats(src),
                                 readSfxBWTSym, params,
                                 numExtHeaders, headerIDs, extHeaderSizes,
                                 extHeaderCallbacks, headerCBData, biFunc,
                                 cwExtBitsPerPos, biVarBits,
                                 cbState, err);
  if (!newSeqIdx)
    MRAEncDelete(alphabet);
  return newSeqIdx;
}

extern EISeq *
createEncIdxSeqGen(Seqpos totalLen, const Str *projectName,
                   MRAEnc *alphabet, const struct seqStats *stats,
                   SeqDataReader seqGenerator,
                   const struct seqBaseParam *params,
                   size_t numExtHeaders, const uint16_t *headerIDs,
                   const uint32_t *extHeaderSizes,
                   headerWriteFunc *extHeaderCallbacks,
                   void **headerCBData,
                   bitInsertFunc biFunc, BitOffset cwExtBitsPerPos,
                   varExtBitsEstimator biVarBits, void *cbState, Error *err)
{
  EISeq *seqIdx = NULL;
  switch (params->encType)
  {
  case BWT_ON_BLOCK_ENC:
    seqIdx = newGenBlockEncIdxSeq(totalLen, projectName, alphabet, stats,
                                  seqGenerator, params, numExtHeaders,
                                  headerIDs, extHeaderSizes, extHeaderCallbacks,
                                  headerCBData, biFunc, cwExtBitsPerPos,
                                  biVarBits, cbState, err);
    break;
  default:
    error_set(err, "Illegal/unknown/unimplemented encoding requested!");
    break;
  }
  return seqIdx;
}

extern struct encIdxSeq *
loadEncIdxSeqForSA(const Suffixarray *sa, Seqpos totalLen,
                   const Str *projectName,
                   enum seqBaseEncoding encType, int features, Error *err)
{
  MRAEnc *alphabet;
  EISeq *seqIdx = NULL;
  assert(sa);
  alphabet = SANewMRAEnc(sa);
  switch (encType)
  {
  case BWT_ON_BLOCK_ENC:
    seqIdx = loadBlockEncIdxSeqGen(alphabet, totalLen, projectName, features,
                                   err);
    break;
  default:
    error_set(err, "Illegal/unknown/unimplemented encoding requested!");
    break;
  }
  if (!seqIdx)
    MRAEncDelete(alphabet);
  return seqIdx;
}

extern EISeq *
loadEncIdxSeq(const Str *projectName,
              enum seqBaseEncoding encType, int features,
              Verboseinfo *verbosity, Error *err)
{
  struct encIdxSeq *newSeqIdx = NULL;
  Suffixarray suffixArray;
  Seqpos len;
  do
  {
    if (streamsuffixarray(&suffixArray, &len,
                          0, projectName, verbosity, err))
      break;
    ++len;
    newSeqIdx = loadEncIdxSeqForSA(&suffixArray, len, projectName,
                                   encType, features, err);
    freesuffixarray(&suffixArray);
  } while (0);
  return newSeqIdx;
}
