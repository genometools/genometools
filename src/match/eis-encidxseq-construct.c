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

#include "core/encseq_metadata.h"
#include "match/eis-blockcomp-construct.h"
#include "match/eis-encidxseq-construct.h"
#include "match/sarr-def.h"
#include "match/esa-map.h"
#include "core/log_api.h"

static EISeq *
gt_createEncIdxSeqFromSASeqSrc(SASeqSrc *src,
                            const char *projectName,
                            const struct seqBaseParam *params,
                            size_t numExtHeaders, const uint16_t *headerIDs,
                            const uint32_t *extHeaderSizes,
                            headerWriteFunc *extHeaderCallbacks,
                            void **headerCBData, bitInsertFunc biFunc,
                            BitOffset cwExtBitsPerPos,
                            varExtBitsEstimator biVarBits, void *cbState,
                            GtError *err);

EISeq *
gt_createEncIdxSeq(const char *projectName,
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
  unsigned long length;
  gt_assert(projectName);
  /* map and interpret index project file */
  if (streamsuffixarray(&suffixArray,
                       SARR_SUFTAB | SARR_BWTTAB, projectName, verbosity, err))
    return NULL;
  length = gt_encseq_total_length(suffixArray.encseq) + 1;
  newSeqIdx = gt_createEncIdxSeqFromSA(&suffixArray, length,
                                      projectName, params,
                                      numExtHeaders, headerIDs,
                                      extHeaderSizes, extHeaderCallbacks,
                                      headerCBData, biFunc, cwExtBitsPerPos,
                                      biVarBits, cbState, err);
  gt_freesuffixarray(&suffixArray);
  return newSeqIdx;
}

EISeq *
gt_createEncIdxSeqFromSA(Suffixarray *sa, unsigned long totalLen,
                      const char *projectName,
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
  gt_initSuffixarrayFileInterface(&sai, totalLen, sa);
  newSeqIdx = gt_createEncIdxSeqFromSAI(
    &sai, projectName, params, numExtHeaders, headerIDs,
    extHeaderSizes, extHeaderCallbacks, headerCBData, biFunc, cwExtBitsPerPos,
    biVarBits, cbState, err);
  gt_destructSuffixarrayFileInterface(&sai);
  return newSeqIdx;
}

EISeq *
gt_createEncIdxSeqFromSAI(SuffixarrayFileInterface *sai,
                       const char *projectName,
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
  return gt_createEncIdxSeqFromSASeqSrc(
    SAI2SASS(sai), projectName, params,
    numExtHeaders, headerIDs, extHeaderSizes,
    extHeaderCallbacks, headerCBData, biFunc,
    cwExtBitsPerPos, biVarBits, cbState, err);
}

EISeq *
gt_createEncIdxSeqFromSfxI(sfxInterface *sfxi,
                        const char *projectName,
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
  return gt_createEncIdxSeqFromSASeqSrc(
    gt_SfxI2SASS(sfxi), projectName, params,
    numExtHeaders, headerIDs, extHeaderSizes,
    extHeaderCallbacks, headerCBData, biFunc,
    cwExtBitsPerPos, biVarBits, cbState, err);
}

static EISeq *
gt_createEncIdxSeqFromSASeqSrc(SASeqSrc *src,
                            const char *projectName,
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
  newSeqIdx = gt_createEncIdxSeqGen(SASSGetLength(src), projectName,
                                 alphabet, SASSGetSeqStats(src),
                                 readSfxBWTSym, params,
                                 numExtHeaders, headerIDs, extHeaderSizes,
                                 extHeaderCallbacks, headerCBData, biFunc,
                                 cwExtBitsPerPos, biVarBits,
                                 cbState, err);
  if (!newSeqIdx)
    gt_MRAEncDelete(alphabet);
  return newSeqIdx;
}

EISeq *
gt_createEncIdxSeqGen(unsigned long totalLen, const char *projectName,
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
    seqIdx = gt_newGenBlockEncIdxSeq(totalLen, projectName, alphabet, stats,
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

struct encIdxSeq *
gt_loadEncIdxSeqForSA(const GtAlphabet *gtalphabet,
                   const char *projectName,
                   enum seqBaseEncoding encType, int features, GtError *err)
{
  MRAEnc *alphabet;
  EISeq *seqIdx = NULL;
  gt_assert(gtalphabet!=NULL);
  alphabet = gt_SANewMRAEnc(gtalphabet);
  switch (encType)
  {
  case BWT_ON_BLOCK_ENC:
    seqIdx = gt_loadBlockEncIdxSeqGen(alphabet, projectName, features,
                                   err);
    break;
  default:
    gt_error_set(err, "Illegal/unknown/unimplemented encoding requested!");
    break;
  }
  gt_assert(seqIdx != NULL || gt_error_is_set(err));
  return seqIdx;
}

EISeq *
gt_loadEncIdxSeq(const char *projectName,
                 enum seqBaseEncoding encType, int features,
                 GT_UNUSED GtLogger *verbosity, GtError *err)
{
  struct encIdxSeq *newSeqIdx = NULL;
  GtEncseq *encseq = NULL;
  GtEncseqLoader *el = NULL;
  do
  {
    el = gt_encseq_loader_new();
    gt_encseq_loader_do_not_require_sds_tab(el);
    gt_encseq_loader_do_not_require_des_tab(el);
    gt_encseq_loader_do_not_require_ssp_tab(el);
    encseq = gt_encseq_loader_load(el, projectName, err);
    gt_encseq_loader_delete(el);
    if (encseq == NULL)
      break;
    newSeqIdx = gt_loadEncIdxSeqForSA(gt_encseq_alphabet(encseq),
                                      projectName,
                                      encType, features, err);
    gt_encseq_delete(encseq);
  } while (0);
  return newSeqIdx;
}
