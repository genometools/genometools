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
#ifndef EIS_SEQDATASRC_H
#define EIS_SEQDATASRC_H
/**
 * \file eis-seqdatasrc.h
 * Abstract interface for sequential data sources.
 */
#include <stdlib.h>
#include <string.h>

#include "libgtcore/error.h"

/* sequential data source api */
typedef void *SeqDataSrc;
/**
 * @param len number of elements (objects) to read
 * @return number of elements read (less then len on end-of-file etc.)
 */
typedef size_t (*seqDataReadFunc)(SeqDataSrc src, void *dest, size_t len,
                                  Error *err);

struct seqDataReader
{
  SeqDataSrc src;
  seqDataReadFunc readData;
};

typedef struct seqDataReader SeqDataReader;

static inline int
SDRIsValid(SeqDataReader sr)
{
  return sr.readData != NULL;
}

static inline size_t
SDRRead(SeqDataReader sr, void *dest, size_t len, Error *err)
{
  return sr.readData(sr.src, dest, len, err);
}

/* while the above describes a data source, this is a corresponding sink */
typedef void *SeqDataDest;
typedef size_t (*seqDataWriteFunc)(SeqDataDest dest, const void *src,
                                   size_t len);

struct seqDataWriter
{
  SeqDataDest dest;
  seqDataWriteFunc writeData;
};

typedef struct seqDataWriter SeqDataWriter;

static inline size_t
SDWWrite(SeqDataWriter sw, const void *src, size_t len)
{
  return sw.writeData(sw.dest, src, len);
}

/* generic data translator api */
union translatorState
{
  void *ref;
  size_t elemSize;
};

typedef union translatorState TranslatorState;
/**
 * @return number of chars! written to dest
 */
typedef size_t (*seqDataTranslateFunc)(void *translator, void *dest,
                                       const void *src, size_t len);

struct seqDataTranslator
{
  TranslatorState state;
  seqDataTranslateFunc translateData;
};

typedef struct seqDataTranslator SeqDataTranslator;

static inline size_t
SDRTranslate(SeqDataTranslator xltor, void *dest, const void *src, size_t len)
{
  if (xltor.translateData)
    return xltor.translateData(xltor.state.ref, dest, src, len);
  /* fall back to zero-translation i.e. verbatim copy */
  memcpy(dest, src, len * xltor.state.elemSize);
  return len * xltor.state.elemSize;
}

#endif
