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

#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "match/dataalign.h"
#include "core/fa.h"
#include "core/filelengthvalues.h"
#include "core/logger.h"
#include "core/minmax.h"
#include "core/seq_iterator_api.h"
#include "core/str.h"
#include "core/str_array.h"
#include "core/types_api.h"
#include "core/unused_api.h"
#include "core/codetype.h"
#include "core/encseq.h"

#include "match/eis-encidxseq.h"
#include "match/eis-sa-common.h"
#include "match/eis-sequencemultiread.h"
#include "match/eis-suffixerator-interface.h"

struct sfxInterface
{
  struct SASeqSrc baseClass;
  GtReadmode readmode;
  unsigned int prefixlength, numofparts, maximumspace;
  const Sfxstrategy *sfxstrategy;
  const GtAlphabet *alpha;
  const GtEncseq *encseq;
  struct seqStats *stats;
  Sfxiterator *sfi;
  Definedunsignedlong rot0Pos;
  /* data relevant to holding portions of the suffix array */
  unsigned long lastGeneratedLen, lastGeneratedStart;
  const GtSuffixsortspace *lastGeneratedSufTabSegment;
};

static SeqDataTranslator
SfxIRequest2XltorFunc(sfxInterface *sfxi,
                      enum sfxDataRequest rtype)
{
  SeqDataTranslator tr = { { NULL }, NULL };
  switch (rtype)
  {
    struct encSeqTrState readState;
    struct saTaggedXltorState *stateStore;
  case SFX_REQUEST_BWTTAB:
    readState.readmode = sfxi->readmode;
    readState.encseq = sfxi->encseq;
    stateStore = gt_addSuffixarrayXltor(&sfxi->baseClass.xltorStates,
                                        rtype, readState);
    tr.state.ref = &stateStore->state;
    tr.translateData = gt_translateSuftab2BWT;
    tr.translateDataSuffixsortspace = gt_translateSuftab2BWTSuffixsortspace;
    break;
  case SFX_REQUEST_SUFTAB:
    tr.state.elemSize = sizeof (unsigned long);
    gt_assert(tr.translateData == NULL);
    break;
  default:
    fprintf(stderr, "error: unimplemented request!\n");
    abort();
    break;
  }
  return tr;
}

static inline sfxInterface *
SASS2SfxI(SASeqSrc *baseClass)
{
  return (sfxInterface *)((char *)baseClass
                          - offsetof(sfxInterface, baseClass));
}

struct SASeqSrc *
gt_SfxI2SASS(sfxInterface *sfxi)
{
  return &sfxi->baseClass;
}

static inline const sfxInterface *
constSASS2SfxI(const SASeqSrc *baseClass)
{
  return (const sfxInterface *)((const char *)baseClass
                                - offsetof(sfxInterface, baseClass));
}

static SeqDataTranslator
SfxIBaseRequest2XltorFunc(SASeqSrc *baseClass,
                          enum sfxDataRequest rtype)
{
  return SfxIRequest2XltorFunc(SASS2SfxI(baseClass), rtype);
}

static Definedunsignedlong
SfxIBaseGetRot0Pos(const SASeqSrc *baseClass)
{
  return gt_SfxIGetRot0Pos(constSASS2SfxI(baseClass));
}

static const struct seqStats *
SfxIBaseGetSeqStats(const SASeqSrc *baseClass)
{
  return gt_SfxIGetSeqStats(constSASS2SfxI(baseClass));
}

static MRAEnc *
SfxIBaseNewMRAEnc(const SASeqSrc *baseClass)
{
  return gt_SfxINewMRAEnc(constSASS2SfxI(baseClass));
}

static void
gt_deleteSfxInterfaceBase(SASeqSrc *baseClass)
{
  gt_deleteSfxInterface(SASS2SfxI(baseClass));
}

static size_t
SfxIGenerate(void *iface,
             void *backlogState,
             move2BacklogFunc move2Backlog,
             void *output,
             unsigned long generateStart,
             size_t len,
             SeqDataTranslator xltor);

sfxInterface *
gt_newSfxInterface(GtReadmode readmode,
                unsigned int prefixlength,
                unsigned int numofparts,
                unsigned long maximumspace,
                const Sfxstrategy *sfxstrategy,
                const GtEncseq *encseq,
                GtTimer *sfxprogress,
                bool withprogressbar,
                unsigned long length,
                GtLogger *verbosity,
                GtError *err)
{
  return gt_newSfxInterfaceWithReaders(readmode,
                                       prefixlength,
                                       numofparts,
                                       maximumspace,
                                       sfxstrategy,
                                       0,
                                       NULL,
                                       NULL,
                                       encseq,
                                       sfxprogress,
                                       withprogressbar,
                                       length,
                                       verbosity,
                                       err);
}

static struct seqStats *
newSeqStatsFromCharDist(const GtEncseq *encseq,
                        const GtAlphabet *alpha, unsigned long len)
{
  struct seqStats *stats = NULL;
  unsigned i, numofchars;
  unsigned long regularSymsSum = 0;
  stats = gt_malloc(offsetAlign(sizeof (*stats), sizeof (unsigned long))
                    + (UINT8_MAX + 1) * sizeof (unsigned long));
  unsigned int numOfSeqs;

  numOfSeqs = gt_encseq_num_of_sequences(encseq);
  stats->sourceAlphaType = sourceUInt8;
  stats->symbolDistributionTable =
    (unsigned long *)((char *)stats + offsetAlign(sizeof (*stats),
                                                  sizeof (unsigned long)));
  memset(stats->symbolDistributionTable,
         0,
         sizeof (unsigned long) * (UINT8_MAX + 1));
  numofchars = gt_alphabet_num_of_chars(alpha);
  for (i = 0; i < numofchars; ++i)
  {
    stats->symbolDistributionTable[i]
      = (unsigned long) gt_encseq_charcount(encseq,(GtUchar) i);
    regularSymsSum += stats->symbolDistributionTable[i];
  }
  stats->symbolDistributionTable[WILDCARD] = len - regularSymsSum - numOfSeqs;
  stats->symbolDistributionTable[SEPARATOR] += numOfSeqs;
  stats->symbolDistributionTable[UNDEFBWTCHAR] += 1;
  return stats;
}

static void
deleteSeqStats(struct seqStats *stats)
{
  gt_free(stats);
}

#define gt_newSfxInterfaceWithReadersErrRet()        \
  do {                                            \
    if (sfxi->stats)                             \
      deleteSeqStats(sfxi->stats);               \
    if (sfxi) gt_free(sfxi);                    \
    sfxi = NULL;                                 \
  } while (0)

sfxInterface *
gt_newSfxInterfaceWithReaders(GtReadmode readmode,
                              unsigned int prefixlength,
                              unsigned int numofparts,
                              unsigned long maximumspace,
                              const Sfxstrategy *sfxstrategy,
                              size_t numReaders,
                              enum sfxDataRequest readerRequests[],
                              SeqDataReader readers[],
                              const GtEncseq *encseq,
                              GtTimer *sfxprogress,
                              bool withprogressbar,
                              unsigned long length,
                              GtLogger *verbosity, GtError *err)
{
  sfxInterface *sfxi = NULL;
  size_t idx;

  gt_error_check(err);
  sfxi = gt_calloc(1, sizeof (*sfxi));
  {
    RandomSeqAccessor origSeqAccess = { gt_SfxIGetOrigSeq, sfxi };
    initSASeqSrc(&sfxi->baseClass, length, SfxIBaseRequest2XltorFunc, NULL,
                 SfxIBaseGetRot0Pos, SfxIBaseGetSeqStats,
                 origSeqAccess, gt_deleteSfxInterfaceBase, SfxIBaseNewMRAEnc,
                 SfxIGenerate, sfxi);
  }
  sfxi->readmode = readmode;
  sfxi->encseq = encseq;
  sfxi->alpha = gt_encseq_alphabet(encseq);
  sfxi->stats = newSeqStatsFromCharDist(encseq,sfxi->alpha, length);
  if (!(sfxi->sfi = gt_Sfxiterator_new(encseq,
                                       readmode,
                                       prefixlength,
                                       numofparts,
                                       maximumspace,
                                       sfxstrategy,
                                       sfxprogress,
                                       withprogressbar,
                                       verbosity,
                                       err)))
  {
    gt_newSfxInterfaceWithReadersErrRet();
  }
  sfxi->rot0Pos.defined = false;

  sfxi->lastGeneratedStart = sfxi->lastGeneratedLen = 0;
  sfxi->lastGeneratedSufTabSegment = NULL;

  for (idx = 0; idx < numReaders; ++idx)
  {
    readers[idx] = gt_SfxIRegisterReader(sfxi, readerRequests[idx]);
    if (!readers[idx].readData)
    {
      gt_newSfxInterfaceWithReadersErrRet();
    }
  }
  return sfxi;
}

const Sfxiterator *gt_SfxInterface2Sfxiterator(const sfxInterface *sfxi)
{
  return sfxi->sfi;
}

void
gt_deleteSfxInterface(sfxInterface *sfxi)
{
  destructSASeqSrc(&sfxi->baseClass);
  (void) gt_Sfxiterator_delete(sfxi->sfi,NULL);
  sfxi->sfi = NULL;
  deleteSeqStats(sfxi->stats);
  gt_free(sfxi);
}

const GtAlphabet *
gt_SfxIGetAlphabet(const sfxInterface *si)
{
  return si->alpha;
}

MRAEnc *
gt_SfxINewMRAEnc(const sfxInterface *si)
{
  MRAEnc *alphabet;
  gt_assert(si);
  alphabet = gt_MRAEncGTAlphaNew(gt_SfxIGetAlphabet(si));
  gt_MRAEncAddSymbolToRange(alphabet, SEPARATOR, 1);
  return alphabet;
}

unsigned long
gt_SfxIGetLength(const sfxInterface *si)
{
  gt_assert(si);
  return si->baseClass.seqLen;
}

const struct seqStats *
gt_SfxIGetSeqStats(const sfxInterface *si)
{
  return si->stats;
}

Definedunsignedlong
gt_SfxIGetRot0Pos(const struct sfxInterface *si)
{
  return si->rot0Pos;
}

const GtEncseq *
gt_SfxIGetEncSeq(const sfxInterface *si)
{
  return si->encseq;
}

GtReadmode
gt_SfxIGetReadmode(const sfxInterface *si)
{
  return si->readmode;
}

SeqDataReader
gt_SfxIRegisterReader(sfxInterface *sfxi, enum sfxDataRequest rtype)
{
  return gt_seqReaderSetRegisterConsumer(
    &sfxi->baseClass.readerSet, rtype, SfxIRequest2XltorFunc(sfxi, rtype));
}

size_t
gt_SfxIGetOrigSeq(const void *state, Symbol *dest, unsigned long pos,
                  size_t len)
{
  const struct sfxInterface *sfxi;
  gt_assert(state);
  sfxi = state;
  return EncSeqGetSubSeq(sfxi->encseq, sfxi->readmode, pos, len, dest);
}

/** writes substring of suffix table to output, stores older data into
 * cache if necessary */
static size_t
SfxIGenerate(void *iface,
             void *backlogState,
             move2BacklogFunc move2Backlog,
             void *output,
             unsigned long generateStart,
             size_t len,
             SeqDataTranslator xltor)
{
  sfxInterface *sfxi = iface;
  size_t elemsLeft = len;
  gt_assert(sfxi && backlogState && move2Backlog && output);
  gt_assert(generateStart + len <= gt_SfxIGetLength(sfxi));
  do
  {
    if (generateStart < sfxi->lastGeneratedStart + sfxi->lastGeneratedLen)
    {
      size_t copyLen = MIN(elemsLeft, sfxi->lastGeneratedStart
                           + sfxi->lastGeneratedLen - generateStart),
        charsWritten =
        SDRTranslateSuffixsortspace(xltor, output,
                                    sfxi->lastGeneratedSufTabSegment,
                                    generateStart - sfxi->lastGeneratedStart,
                                    copyLen);
      generateStart += copyLen;
      elemsLeft -= copyLen;
      output = (char *)output + charsWritten;
    }
    /* 1. read next region of sequence by calling nextSfxIterator */
    if (elemsLeft)
    {

      move2Backlog(backlogState, sfxi->lastGeneratedSufTabSegment,
                   sfxi->lastGeneratedStart, sfxi->lastGeneratedLen);
      sfxi->lastGeneratedStart += sfxi->lastGeneratedLen;
      sfxi->lastGeneratedSufTabSegment
        = gt_Sfxiterator_next(&sfxi->lastGeneratedLen,
                              NULL,
                              sfxi->sfi);
      if (sfxi->lastGeneratedSufTabSegment != NULL)
      {
        /* size_t because the current approach cannot generate more
         * than memory will hold anyway */
        size_t pos, lastGeneratedLen = sfxi->lastGeneratedLen;

        if (!sfxi->rot0Pos.defined)
        {
          for (pos=0; pos < lastGeneratedLen; pos++)
          {
            if (gt_suffixsortspace_getdirect(sfxi->lastGeneratedSufTabSegment,
                                             pos) == 0)
            {
              sfxi->rot0Pos.defined = true;
              sfxi->rot0Pos.valueunsignedlong = sfxi->lastGeneratedStart + pos;
              break;
            }
          }
        }
        /* uncomment this to reenable synchronous writing of tables */
/*if (sfxi->lastGeneratedSufTabSegment == NULL */
/*    || suftab2file(&sfxi->outfileinfo, sfxi->lastGeneratedSufTabSegment, */
/*                   sfxi->so.readmode, sfxi->lastGeneratedLen, err) != 0) */
/*       break; */
      }
      else
        break;
    }
    /* 5. if positions in region don't suffice go back to step 3. */
  } while (elemsLeft);
  return len - elemsLeft;
}
