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

#include <stdlib.h>

#include "core/chardef.h"
#include "core/error.h"
#include "core/ma_api.h"
#include "core/unused_api.h"
#include "match/sarr-def.h"
#include "match/eis-mrangealphabet.h"
#include "match/eis-sa-common.h"
#include "match/eis-seqdatasrc.h"
#include "match/eis-suffixarray-interface.h"
#include "match/eis-suffixarray-interface-priv.h"

static void
gt_deleteSuffixarrayFileInterfaceBase(SASeqSrc *baseClass)
{
  gt_deleteSuffixarrayFileInterface(SASS2SAI(baseClass));
}

static struct seqDataReader
SAIBaseMakeReader(SASeqSrc *baseClass, enum sfxDataRequest rtype)
{
  return gt_SAIMakeReader(SASS2SAI(baseClass), rtype);
}

static Definedunsignedlong
SAIBaseGetRot0Pos(const SASeqSrc *baseClass)
{
  return gt_SAIGetRot0Pos(constSASS2SAI(baseClass));
}

MRAEnc *
gt_SAIBaseNewMRAEnc(const SASeqSrc *baseClass)
{
  return SAINewMRAEnc(constSASS2SAI(baseClass));
}

static size_t
SAIGenerate(void *generatorState, void *backlogState,
            move2BacklogFunc move2Backlog, void *output,
            GtUword generateStart, size_t len,
            SeqDataTranslator xltor);

void
gt_initSuffixarrayFileInterface(SuffixarrayFileInterface *sai,
                             GtUword seqLen, Suffixarray *sa)
{
  {
    RandomSeqAccessor origSeqAccess = { gt_SAIGetOrigSeq, sai };
    initSASeqSrc(&sai->baseClass, seqLen, NULL, SAIBaseMakeReader,
                 SAIBaseGetRot0Pos, NULL,
                 origSeqAccess, gt_deleteSuffixarrayFileInterfaceBase,
                 gt_SAIBaseNewMRAEnc,
                 SAIGenerate, sai);
  }
  sai->sa = sa;
  sai->numBWTFileReaders = 0;
  gt_initSATaggedXltorStateList(&sai->xltorStates);
}

SuffixarrayFileInterface *
gt_newSuffixarrayFileInterface(Suffixarray *sa, GtUword seqLen)
{
  SuffixarrayFileInterface *sai = gt_malloc(sizeof (*sai));
  gt_initSuffixarrayFileInterface(sai, seqLen, sa);
  return sai;
}

void
gt_destructSuffixarrayFileInterface(SuffixarrayFileInterface *sai)
{
  destructSASeqSrc(&sai->baseClass);
  gt_destructSATaggedXltorStateList(&sai->xltorStates);
}

void
gt_deleteSuffixarrayFileInterface(SuffixarrayFileInterface *sai)
{
  gt_destructSuffixarrayFileInterface(sai);
  gt_free(sai);
}

static size_t
SAIReadBWT(void *state, Symbol *dest, size_t len, GtError *err);

struct seqDataReader
gt_SAIMakeBWTReader(SuffixarrayFileInterface *sai)
{
  struct seqDataReader reader = { NULL, NULL};
  if (!sai->sa->bwttabstream.fp || sai->numBWTFileReaders > 0)
  {
    if (sai->sa->encseq)
    {
      struct encSeqTrState bwtReadState = {
        .readmode = sai->sa->readmode,
        .encseq = sai->sa->encseq
      };
      struct saTaggedXltorState *stateStore
        = gt_addSuffixarrayXltor(&sai->xltorStates,
                              SFX_REQUEST_BWTTAB, bwtReadState);
      struct seqDataTranslator xltor = {
        { .ref = &stateStore->state },
        gt_translateSuftab2BWT,
        gt_translateSuftab2BWTSuffixsortspace
      };

      reader = gt_seqReaderSetRegisterConsumer(&sai->baseClass.readerSet,
                                            SFX_REQUEST_BWTTAB, xltor);
    }
    else
    {
      fputs("error: bwt data not available for given project.\n", stderr);
    }
  }
  else
  {
    /* a .bwt file is available for reading */
    reader.readData = (seqDataReadFunc)SAIReadBWT;
    reader.src = sai;
    ++sai->numBWTFileReaders;
  }
  return reader;
}

struct seqDataReader
gt_SAIMakeSufTabReader(SuffixarrayFileInterface *sai)
{
  struct seqDataReader reader = { NULL, NULL};
  if (sai->sa->suftabstream_GtUlong.fp)
  {
    struct seqDataTranslator xltor = {
      { .elemSize = sizeof (GtUword) }, NULL, NULL,
    };
    reader = gt_seqReaderSetRegisterConsumer(&sai->baseClass.readerSet,
                                             SFX_REQUEST_SUFTAB, xltor);
  }
  else
  {
    fputs("error: suffix array data not available for given project.\n",
          stderr);
  }
  return reader;
}

struct seqDataReader
gt_SAIMakeReader(SuffixarrayFileInterface *sai, enum sfxDataRequest rtype)
{
  struct seqDataReader reader = { NULL, NULL};
  switch (rtype)
  {
  case SFX_REQUEST_SUFTAB:
    reader = gt_SAIMakeSufTabReader(sai);
    break;
  case SFX_REQUEST_BWTTAB:
    reader = gt_SAIMakeBWTReader(sai);
    break;
  default:
    fprintf(stderr, "error: unimplemented request: %d, %s: %d!\n", rtype,
            __FILE__, __LINE__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  return reader;
}

/**
 * @brief Read given length of symbols from the BWT, starting after last
 * position read.
 * @param state reference of a SuffixarrayFileInterface
 * @param dest write symbols here
 * @param len length of string to read
 * @return actual number of symbols read
 */
static size_t
SAIReadBWT(void *state, GtUchar *dest, size_t len, GT_UNUSED GtError *err)
{
  SuffixarrayFileInterface *sai = state;
  gt_assert(state);
  return fread(dest, sizeof (GtUchar), len, sai->sa->bwttabstream.fp);
}

size_t
gt_SAIGetOrigSeq(const void *state, Symbol *dest, GtUword pos, size_t len)
{
  const SuffixarrayFileInterface *sai;
  gt_assert(state);
  sai = state;
  return EncSeqGetSubSeq(sai->sa->encseq, sai->sa->readmode, pos, len, dest);
}

Definedunsignedlong
gt_SAIGetRot0Pos(const void *state)
{
  const SuffixarrayFileInterface *sai = state;
  gt_assert(sai);
  return sai->sa->longest;
}

MRAEnc *
gt_SANewMRAEnc(const GtAlphabet *gtalphabet)
{
  MRAEnc *alphabet;
  gt_assert(gtalphabet != NULL);
  alphabet = gt_MRAEncGTAlphaNew(gtalphabet);
  gt_MRAEncAddSymbolToRange(alphabet, SEPARATOR, 1);
  return alphabet;
}

static size_t
SAIGenerate(void *generatorState, void *backlogState,
            move2BacklogFunc move2Backlog, void *output,
            GtUword generateStart, size_t len,
            SeqDataTranslator xltor)
{
  size_t idx;
  SuffixarrayFileInterface *sai = generatorState;
  Suffixarray *sa;
  GtUword buf[len];

  gt_assert(sai);
  sa = sai->sa;
  for (idx = 0; idx < len; ++idx)
  {
    if (gt_readnextfromstream_GtUlong(buf + idx,&sa->suftabstream_GtUlong) != 1)
    {
      break;
    }
  }
  move2Backlog(backlogState, buf, generateStart, idx);
  SDRTranslate(xltor, output, buf, idx);
  return idx;
}
