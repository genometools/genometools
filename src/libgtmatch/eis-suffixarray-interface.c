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

#include "libgtcore/chardef.h"
#include "libgtcore/error.h"
#include "libgtcore/ma.h"
#include "libgtcore/unused.h"
#include "libgtmatch/sarr-def.h"
#include "libgtmatch/eis-mrangealphabet.h"
#include "libgtmatch/eis-sa-common.h"
#include "libgtmatch/eis-seqdatasrc.h"
#include "libgtmatch/eis-suffixarray-interface.h"
#include "libgtmatch/eis-suffixarray-interface-priv.h"

static void
deleteSuffixarrayFileInterfaceBase(SASeqSrc *baseClass)
{
  deleteSuffixarrayFileInterface(SASS2SAI(baseClass));
}

static struct seqDataReader
SAIBaseMakeReader(SASeqSrc *baseClass, enum sfxDataRequest rtype)
{
  return SAIMakeReader(SASS2SAI(baseClass), rtype);
}

static DefinedSeqpos
SAIBaseGetRot0Pos(const SASeqSrc *baseClass)
{
  return SAIGetRot0Pos(constSASS2SAI(baseClass));
}

static inline size_t
SAIBaseGetOrigSeq(const SASeqSrc *baseClass, Symbol *dest, Seqpos pos,
                  size_t len)
{
  return SAIGetOrigSeq(constSASS2SAI(baseClass), dest, pos, len);
}

extern MRAEnc *
SAIBaseNewMRAEnc(const SASeqSrc *baseClass)
{
  return SAINewMRAEnc(constSASS2SAI(baseClass));
}

static size_t
SAIGenerate(void *generatorState, void *backlogState,
            move2BacklogFunc move2Backlog, void *output,
            Seqpos generateStart, size_t len,
            SeqDataTranslator xltor, Error *err);

extern void
initSuffixarrayFileInterface(SuffixarrayFileInterface *sai, Seqpos seqLen,
                             Suffixarray *sa)
{
  {
    RandomSeqAccessor origSeqAccess = { SAIGetOrigSeq, sai };
    initSASeqSrc(&sai->baseClass, seqLen, NULL, SAIBaseMakeReader,
                 SAIBaseGetRot0Pos, NULL,
                 origSeqAccess, deleteSuffixarrayFileInterfaceBase,
                 SAIBaseNewMRAEnc,
                 SAIGenerate, sai);
  }
  sai->sa = sa;
  sai->numBWTFileReaders = 0;
  initSATaggedXltorStateList(&sai->xltorStates);
}

extern SuffixarrayFileInterface *
newSuffixarrayFileInterface(Suffixarray *sa, Seqpos seqLen)
{
  SuffixarrayFileInterface *sai = ma_malloc(sizeof (*sai));
  initSuffixarrayFileInterface(sai, seqLen, sa);
  return sai;
}

extern void
destructSuffixarrayFileInterface(SuffixarrayFileInterface *sai)
{
  destructSASeqSrc(&sai->baseClass);
  destructSATaggedXltorStateList(&sai->xltorStates);
}

extern void
deleteSuffixarrayFileInterface(SuffixarrayFileInterface *sai)
{
  destructSuffixarrayFileInterface(sai);
  ma_free(sai);
}

static size_t
SAIReadBWT(void *state, Symbol *dest, size_t len, Error *err);

extern struct seqDataReader
SAIMakeBWTReader(SuffixarrayFileInterface *sai)
{
  struct seqDataReader reader = { NULL, NULL};
  if (!sai->sa->bwttabstream.fp || sai->numBWTFileReaders > 0)
  {
    if (sai->sa->encseq)
    {
      union saXltorState bwtReadState = {
        .encSeqTr.readmode = sai->sa->readmode,
        .encSeqTr.encseq = sai->sa->encseq
      };
      struct saTaggedXltorState *stateStore
        = addSuffixarrayXltor(&sai->xltorStates,
                              SFX_REQUEST_BWTTAB, bwtReadState);
      struct seqDataTranslator xltor = {
        { .ref = &stateStore->state.encSeqTr },
        (seqDataTranslateFunc)translateSuftab2BWT
      };

      reader = seqReaderSetRegisterConsumer(&sai->baseClass.readerSet,
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

extern struct seqDataReader
SAIMakeSufTabReader(SuffixarrayFileInterface *sai)
{
  struct seqDataReader reader = { NULL, NULL};
  if (sai->sa->suftabstream.fp)
  {
    struct seqDataTranslator xltor = {
      { .elemSize = sizeof (Seqpos) }, NULL
    };
    reader = seqReaderSetRegisterConsumer(&sai->baseClass.readerSet,
                                          SFX_REQUEST_SUFTAB, xltor);
  }
  else
  {
    fputs("error: suffix array data not available for given project.\n",
          stderr);
  }
  return reader;
}

extern struct seqDataReader
SAIMakeLCPTabReader(SuffixarrayFileInterface *sai)
{
  struct seqDataReader reader = { NULL, NULL};
  if (sai->sa->suftabstream.fp)
  {
    union saXltorState lcpReadState = {
      .lcpState.readmode = sai->sa->readmode,
      .lcpState.encseq = sai->sa->encseq,
      .lcpState.lastSufIdx = -1,
    };
    struct saTaggedXltorState *stateStore
      = addSuffixarrayXltor(&sai->xltorStates,
                            SFX_REQUEST_LCPTAB, lcpReadState);
    struct seqDataTranslator xltor = {
      { .ref = &stateStore->state.lcpState },
      (seqDataTranslateFunc)translateSuftab2BWT
    };
    reader = seqReaderSetRegisterConsumer(&sai->baseClass.readerSet,
                                          SFX_REQUEST_LCPTAB, xltor);
  }
  else
  {
    fputs("error: suffix array data not available for given project.\n",
          stderr);
  }
  return reader;
}

extern struct seqDataReader
SAIMakeReader(SuffixarrayFileInterface *sai, enum sfxDataRequest rtype)
{
  struct seqDataReader reader = { NULL, NULL};
  switch (rtype)
  {
  case SFX_REQUEST_SUFTAB:
    reader = SAIMakeSufTabReader(sai);
    break;
  case SFX_REQUEST_BWTTAB:
    reader = SAIMakeBWTReader(sai);
    break;
  case SFX_REQUEST_LCPTAB:
    reader = SAIMakeLCPTabReader(sai);
    break;
  default:
    fprintf(stderr, "error: unimplemented request: %d, %s: %d!\n", rtype,
            __FILE__, __LINE__);
    exit(EXIT_FAILURE);
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
SAIReadBWT(void *state, Uchar *dest, size_t len, UNUSED Error *err)
{
  SuffixarrayFileInterface *sai = state;
  assert(state);
  return fread(dest, sizeof (Uchar), len, sai->sa->bwttabstream.fp);
}

DECLAREREADFUNCTION(Seqpos)

extern size_t
SAIGetOrigSeq(const void *state, Symbol *dest, Seqpos pos, size_t len)
{
  const SuffixarrayFileInterface *sai;
  assert(state);
  sai = state;
  return EncSeqGetSubSeq(sai->sa->encseq, sai->sa->readmode, pos, len, dest);
}

extern DefinedSeqpos
SAIGetRot0Pos(const void *state)
{
  const SuffixarrayFileInterface *sai = state;
  assert(sai);
  return sai->sa->longest;
}

extern MRAEnc *
SANewMRAEnc(const Suffixarray *sa)
{
  MRAEnc *alphabet;
  assert(sa);
  alphabet = MRAEncGTAlphaNew(sa->alpha);
  MRAEncAddSymbolToRange(alphabet, SEPARATOR, 1);
  return alphabet;
}

static size_t
SAIGenerate(void *generatorState, void *backlogState,
            move2BacklogFunc move2Backlog, void *output,
            Seqpos generateStart, size_t len,
            SeqDataTranslator xltor, Error *err)
{
  size_t i;
  SuffixarrayFileInterface *sai = generatorState;
  Suffixarray *sa;
  Seqpos buf[len];
  assert(sai);
  sa = sai->sa;
  for (i = 0; i < len; ++i)
    if (readnextSeqposfromstream(buf + i, &sa->suftabstream, err) != 1)
      break;
  move2Backlog(backlogState, buf, generateStart, i);
  SDRTranslate(xltor, output, buf, i);
  return i;
}
