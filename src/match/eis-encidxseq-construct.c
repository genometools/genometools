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

#include "match/eis-blockcomp-construct.h"
#include "match/eis-encidxseq-construct.h"
#include "match/sarr-def.h"
#include "match/esa-map.h"

static EISeq *
createEncIdxSeqFromSASeqSrc(SASeqSrc *src,
                            const GtStr *projectName,
                            const struct seqBaseParam *params,
                            size_t numExtHeaders, const uint16_t *headerIDs,
                            const uint32_t *extHeaderSizes,
                            headerWriteFunc *extHeaderCallbacks,
                            void **headerCBData, bitInsertFunc biFunc,
                            BitOffset cwExtBitsPerPos,
                            varExtBitsEstimator biVarBits, void *cbState,
                            GtError *err);

extern EISeq *
createEncIdxSeq(const GtStr *projectName,
                const struct seqBaseParam *params,
                size_t numExtHeaders, const uint16_t *headerIDs,
                const uint32_t *extHeaderSizes,
                headerWriteFunc *extHeaderCallbacks, void **headerCBData,
                bitInsertFunc biFunc, BitOffset cwExtBitsPerPos,
                varExtBitsEstimator biVarBits, void *cbState,
                GtLogger *verbosity, GtError *err)
{
  Suffixarray suffixArray;
  struct encIdxSeq *newSeqIdx;
  Seqpos length;
  gt_assert(projectName);
  /* map and interpret index project file */
  if (streamsuffixarray(&suffixArray,
                       SARR_SUFTAB | SARR_BWTTAB, projectName, verbosity, err))
    return NULL;
  length = gt_encodedsequence_total_length(suffixArray.encseq) + 1;
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
                      const GtStr *projectName,
                      const struct seqBaseParam *params,
                      size_t numExtHeaders, const uint16_t *headerIDs,
                      const uint32_t *extHeaderSizes,
                      headerWriteFunc *extHeaderCallbacks,
                      void **headerCBData,
                      bitInsertFunc biFunc, BitOffset cwExtBitsPerPos,
                      varExtBitsEstimator biVarBits, void *cbState,
                      GtError *err)
{
  struct encIdxSeq *newSeqIdx;
  SuffixarrayFileInterface sai;
  gt_assert(sa && projectName && err);
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
                       const GtStr *projectName,
                       const struct seqBaseParam *params,
                       size_t numExtHeaders, const uint16_t *headerIDs,
                       const uint32_t *extHeaderSizes,
                       headerWriteFunc *extHeaderCallbacks,
                       void **headerCBData,
                       bitInsertFunc biFunc, BitOffset cwExtBitsPerPos,
                       varExtBitsEstimator biVarBits, void *cbState,
                       GtError *err)
{
  gt_assert(sai && projectName && err);
  return createEncIdxSeqFromSASeqSrc(
    SAI2SASS(sai), projectName, params,
    numExtHeaders, headerIDs, extHeaderSizes,
    extHeaderCallbacks, headerCBData, biFunc,
    cwExtBitsPerPos, biVarBits, cbState, err);
}

extern EISeq *
createEncIdxSeqFromSfxI(sfxInterface *sfxi,
                        const GtStr *projectName,
                        const struct seqBaseParam *params,
                        size_t numExtHeaders, const uint16_t *headerIDs,
                        const uint32_t *extHeaderSizes,
                        headerWriteFunc *extHeaderCallbacks,
                        void **headerCBData,
                        bitInsertFunc biFunc, BitOffset cwExtBitsPerPos,
                        varExtBitsEstimator biVarBits, void *cbState,
                        GtError *err)
{
  gt_assert(sfxi && projectName && err);
  return createEncIdxSeqFromSASeqSrc(
    SfxI2SASS(sfxi), projectName, params,
    numExtHeaders, headerIDs, extHeaderSizes,
    extHeaderCallbacks, headerCBData, biFunc,
    cwExtBitsPerPos, biVarBits, cbState, err);
}

static EISeq *
createEncIdxSeqFromSASeqSrc(SASeqSrc *src,
                            const GtStr *projectName,
                            const struct seqBaseParam *params,
                            size_t numExtHeaders, const uint16_t *headerIDs,
                            const uint32_t *extHeaderSizes,
                            headerWriteFunc *extHeaderCallbacks,
                            void **headerCBData, bitInsertFunc biFunc,
                            BitOffset cwExtBitsPerPos,
                            varExtBitsEstimator biVarBits, void *cbState,
                            GtError *err)
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
createEncIdxSeqGen(Seqpos totalLen, const GtStr *projectName,
                   MRAEnc *alphabet, const struct seqStats *stats,
                   SeqDataReader seqGenerator,
                   const struct seqBaseParam *params,
                   size_t numExtHeaders, const uint16_t *headerIDs,
                   const uint32_t *extHeaderSizes,
                   headerWriteFunc *extHeaderCallbacks,
                   void **headerCBData,
                   bitInsertFunc biFunc, BitOffset cwExtBitsPerPos,
                   varExtBitsEstimator biVarBits, void *cbState, GtError *err)
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
    gt_error_set(err, "Illegal/unknown/unimplemented encoding requested!");
    break;
  }
  return seqIdx;
}

extern struct encIdxSeq *
loadEncIdxSeqForSA(const Suffixarray *sa, Seqpos totalLen,
                   const GtStr *projectName,
                   enum seqBaseEncoding encType, int features, GtError *err)
{
  MRAEnc *alphabet;
  EISeq *seqIdx = NULL;
  gt_assert(sa);
  alphabet = SANewMRAEnc(sa);
  switch (encType)
  {
  case BWT_ON_BLOCK_ENC:
    seqIdx = loadBlockEncIdxSeqGen(alphabet, totalLen, projectName, features,
                                   err);
    break;
  default:
    gt_error_set(err, "Illegal/unknown/unimplemented encoding requested!");
    break;
  }
  if (!seqIdx)
    MRAEncDelete(alphabet);
  return seqIdx;
}

extern EISeq *
loadEncIdxSeq(const GtStr *projectName,
              enum seqBaseEncoding encType, int features,
              GtLogger *verbosity, GtError *err)
{
  struct encIdxSeq *newSeqIdx = NULL;
  Suffixarray suffixArray;
  Seqpos len;
  do
  {
    if (streamsuffixarray(&suffixArray, 0, projectName, verbosity, err))
      break;
    len = gt_encodedsequence_total_length(suffixArray.encseq) + 1;
    newSeqIdx = loadEncIdxSeqForSA(&suffixArray, len, projectName,
                                   encType, features, err);
    freesuffixarray(&suffixArray);
  } while (0);
  return newSeqIdx;
}
