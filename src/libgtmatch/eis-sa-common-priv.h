/*
  Copyright (c) 2008 Thomas Jahns <Thomas.Jahns@gmx.net>

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
#ifndef EIS_SA_COMMON_PRIV_H
#define EIS_SA_COMMON_PRIV_H

#include "eis-random-seqaccess.h"
#include "eis-sa-common.h"
#include "eis-sequencemultiread.h"

typedef SeqDataReader (*createReaderFunc)(SASeqSrc *src,
                                          enum sfxDataRequest request);

typedef SeqDataTranslator (*createTranslatorFunc)(SASeqSrc *src,
                                                  enum sfxDataRequest request);

typedef DefinedSeqpos (*getRot0PosFunc)(const SASeqSrc *src);

typedef const struct seqStats *(*getSeqStatsFunc)(const SASeqSrc *src);

typedef void (*deleteSASeqSrcFunc)(SASeqSrc *src);

typedef MRAEnc *(*newMRAEncFunc)(const SASeqSrc *src);

struct SASeqSrc
{
  Seqpos seqLen;
  createReaderFunc createReader;
  createTranslatorFunc createTranslator;
  getRot0PosFunc getRot0Pos;
  getSeqStatsFunc getSeqStats;
  RandomSeqAccessor origSequenceAccess;
  deleteSASeqSrcFunc deleteSASS;
  newMRAEncFunc newMRAEnc;
  MRAEnc *alphabet;
  struct seqReaderSet readerSet;
  struct saTaggedXltorStateList xltorStates;
};

static inline void
initSASeqSrc(SASeqSrc *src, Seqpos seqLen,
             createTranslatorFunc createTranslator,
             createReaderFunc createReader,
             getRot0PosFunc getRot0Pos, getSeqStatsFunc getSeqStats,
             RandomSeqAccessor origSeqAccess, deleteSASeqSrcFunc deleteSASS,
             newMRAEncFunc newSeqMRAEnc,
             generatorFunc generator, void *generatorState);

static inline void
destructSASeqSrc(SASeqSrc *src);

extern SeqDataReader
SASSGenericCreateReader(SASeqSrc *src, enum sfxDataRequest request);

#endif
