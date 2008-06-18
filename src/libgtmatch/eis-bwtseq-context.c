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
#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libgtcore/bitpackstring.h"
#include "libgtcore/dataalign.h"
#include "libgtcore/fa.h"
#include "libgtcore/str.h"
#include "libgtmatch/seqpos-def.h"
#include "libgtmatch/eis-bitpackseqpos.h"
#include "libgtmatch/eis-bwtseq.h"
#include "libgtmatch/eis-bwtseq-context.h"
#include "libgtmatch/eis-bwtseq-context-priv.h"
#include "libgtmatch/eis-seqdatasrc.h"

struct BWTSeqContextRetrieverFactory
{
  Seqpos seqLen, currentSfxPos, moduloMask;
  bool constructionComplete;
  unsigned short mapIntervalLog2;
  FILE *mapTableDiskBackingStore;
  Str *mapTableDBSPath;
};

#define die(msg) die_func(msg, __FILE__, __LINE__)

static inline void
die_func(const char *errMsg, const char *file, int line)
{
  fprintf(stderr, "error: %s: %s at %s, line %d", errMsg, strerror(errno),
          file, line);
  exit(EXIT_FAILURE);
}

enum
{
  BLOCK_IO_SIZE = BUFSIZ / sizeof (Seqpos),
};

static void
initBWTSeqContextRetrieverFactory(BWTSeqContextRetrieverFactory *newFactory,
                                  Seqpos seqLen,
                                  short mapIntervalLog2)
{
  FILE *fp;
  assert(ctxMapILogIsValid(seqLen, mapIntervalLog2));
  if (mapIntervalLog2 == CTX_MAP_ILOG_AUTOSIZE)
  {
    mapIntervalLog2 = requiredUIntBits(requiredSeqposBits(seqLen));
  }
  newFactory->seqLen = seqLen;
  newFactory->currentSfxPos = 0;
  newFactory->moduloMask = (1 << mapIntervalLog2) - 1;
  newFactory->mapIntervalLog2 = mapIntervalLog2;
  newFactory->constructionComplete = false;
  newFactory->mapTableDBSPath = str_new();
  fp = newFactory->mapTableDiskBackingStore
    = fa_xtmpfp_generic(newFactory->mapTableDBSPath,
                        TMPFP_AUTOREMOVE | TMPFP_OPENBINARY);
  {
    off_t backingStoreSize = numMapEntries(seqLen, mapIntervalLog2), i;
    Seqpos buf[BLOCK_IO_SIZE];
    memset(buf, 0, sizeof (buf));
    for (i = BLOCK_IO_SIZE; i < backingStoreSize; i += BLOCK_IO_SIZE)
    {
      if (fwrite(buf, sizeof (buf[0]), BLOCK_IO_SIZE, fp)
          != BLOCK_IO_SIZE)
        die("short write on backing store creation");
    }
    {
      unsigned leftOver = backingStoreSize % BLOCK_IO_SIZE;
      if (leftOver &&
          fwrite(buf, sizeof (buf[0]), leftOver, fp) != leftOver)
        die("short write on backing store creation");
    }
  }
}

extern BWTSeqContextRetrieverFactory *
newBWTSeqContextRetrieverFactory(Seqpos seqLen, short mapIntervalLog2)
{
  BWTSeqContextRetrieverFactory *newFactory;
  newFactory = ma_malloc(sizeof (*newFactory));
  initBWTSeqContextRetrieverFactory(newFactory, seqLen, mapIntervalLog2);
  return newFactory;
}

static void
destructBWTSeqContextRetrieverFactory(BWTSeqContextRetrieverFactory *factory)
{
  fa_xfclose(factory->mapTableDiskBackingStore);
  str_delete(factory->mapTableDBSPath);
}

extern void
deleteBWTSeqContextRetrieverFactory(BWTSeqContextRetrieverFactory *factory)
{
  destructBWTSeqContextRetrieverFactory(factory);
  ma_free(factory);
}

static inline void
addMapVal(BWTSeqContextRetrieverFactory *factory,
          Seqpos currentSfxPos, Seqpos origPos)
{
  off_t mapPos = ((off_t)origPos >> factory->mapIntervalLog2)
    * sizeof (Seqpos);
  Seqpos mapVal = currentSfxPos;
  FILE *fp = factory->mapTableDiskBackingStore;
  if (fseeko(fp, mapPos, SEEK_SET) == -1)
    die("failed to seek in backing store");
  if (fwrite(&mapVal, sizeof (mapVal), 1, fp) != 1)
    die("failed when writing to backing store");
}

/**
 * @return number of not processed suffix indices, i.e. 0 on success
 */
extern Seqpos
BWTSCRFReadAdvance(BWTSeqContextRetrieverFactory *factory, Seqpos chunkSize,
                   SeqDataReader readSfxIdx, Error *err)
{
  Seqpos buf[BLOCK_IO_SIZE], sfxIdxLeft = chunkSize;
  assert(factory);
  while (sfxIdxLeft)
  {
    Seqpos len = MIN(BLOCK_IO_SIZE, sfxIdxLeft);
    if (SDRRead(readSfxIdx, buf, len, err)
        != len)
    {
      fputs("error: short read when building context retriever!\n", stderr);
      abort();
    }
    BWTSCRFMapAdvance(factory, buf, len);
    sfxIdxLeft -= len;
  }
  return chunkSize - sfxIdxLeft;
}

extern size_t
BWTSCRFMapAdvance(BWTSeqContextRetrieverFactory *factory, const Seqpos *src,
                  size_t len)
{
  Seqpos currentSfxPos;
  assert(factory);
  currentSfxPos = factory->currentSfxPos;
  {
    size_t i;
    Seqpos mask = factory->moduloMask,
      seqLen = factory->seqLen;
    for (i = 0; i < len; ++i)
    {
      if (!(((src[i] + seqLen - 1)%seqLen) & mask))
        addMapVal(factory, currentSfxPos + i, ((src[i] + seqLen - 1)%seqLen));
    }
  }
  factory->currentSfxPos = currentSfxPos + len;
  return len;
}

extern bool
BWTSCRFFinished(const BWTSeqContextRetrieverFactory *factory)
{
  return factory->currentSfxPos == factory->seqLen;
}

static inline void
readBS2Map(BWTSeqContextRetrieverFactory *factory,
           BWTSeqContextRetriever *newBWTSeqCR);

static inline bool
BWTSeqCRMapOpen(unsigned short mapIntervalLog2, unsigned short bitsPerSeqpos,
                Seqpos seqLen, const Str *projectName, bool createMapFile,
                BWTSeqContextRetriever *newBWTSeqCR);

extern BWTSeqContextRetriever *
BWTSCRFGet(BWTSeqContextRetrieverFactory *factory, const BWTSeq *bwtSeq,
           const Str *projectName)
{
  unsigned short bitsPerSeqpos, mapIntervalLog2;
  BWTSeqContextRetriever *newBWTSeqCR;
  assert(factory && projectName);
  bitsPerSeqpos = requiredSeqposBits(factory->seqLen - 1);
  newBWTSeqCR = ma_malloc(sizeof (*newBWTSeqCR));
  newBWTSeqCR->mapIntervalLog2 = mapIntervalLog2 = factory->mapIntervalLog2;
  newBWTSeqCR->bitsPerSeqpos = bitsPerSeqpos;
  newBWTSeqCR->bwtSeq = bwtSeq;
  newBWTSeqCR->mapInterval = 1 << factory->mapIntervalLog2;
  newBWTSeqCR->mapMask = newBWTSeqCR->mapInterval - 1;
  if (!BWTSeqCRMapOpen(mapIntervalLog2, bitsPerSeqpos, factory->seqLen,
                       projectName, true, newBWTSeqCR))
  {
    ma_free(newBWTSeqCR);
    return NULL;
  }
  readBS2Map(factory, newBWTSeqCR);
  return newBWTSeqCR;
}

static inline void
readBlock2Buf(FILE *fp, Seqpos *buf, size_t len, BitString bitstring,
              BitOffset offset, unsigned short bitsPerSeqpos)
{
  if (fread(buf, sizeof (buf[0]), len, fp) != len)
    die("short read when reading backing store");
  bsStoreUniformSeqposArray(bitstring, offset, bitsPerSeqpos, len, buf);
}

static inline void
readBS2Map(BWTSeqContextRetrieverFactory *factory,
           BWTSeqContextRetriever *newBWTSeqCR)
{
  FILE *fp = factory->mapTableDiskBackingStore;
  Seqpos buf[BLOCK_IO_SIZE];
  off_t i, numFullBlocks, lastBlockLen;
  BitString revMap =  newBWTSeqCR->revMap;
  BitOffset storePos = 0;
  unsigned bitsPerSeqpos = newBWTSeqCR->bitsPerSeqpos;
  {
    off_t numEntries  = numMapEntries(factory->seqLen,
                                      factory->mapIntervalLog2);
    numFullBlocks = numEntries / BLOCK_IO_SIZE;
    lastBlockLen = numEntries % BLOCK_IO_SIZE;
  }
  if (fseeko(fp, 0, SEEK_SET) == -1)
    die("failed seek in backing store");
  for (i = 0; i < numFullBlocks; ++i )
  {
    readBlock2Buf(fp, buf, BLOCK_IO_SIZE, revMap, storePos, bitsPerSeqpos);
    storePos += bitsPerSeqpos * BLOCK_IO_SIZE;
  }
  readBlock2Buf(fp, buf, lastBlockLen, revMap, storePos, bitsPerSeqpos);
}

enum {
  HEADER_ENTRY_BITS = 16,
};

static inline bool
BWTSeqCRMapOpen(unsigned short mapIntervalLog2, unsigned short bitsPerSeqpos,
                Seqpos seqLen, const Str *projectName, bool createMapFile,
                BWTSeqContextRetriever *newBWTSeqCR)
{
  FILE *mapFile = NULL;
  BitString mapMap = NULL;
  Str *mapName = NULL;
  assert(projectName);
  do {
    size_t headerBitElems = bitElemsAllocSize(2 * HEADER_ENTRY_BITS),
      headerSize = headerBitElems * sizeof (BitElem),
      mapSize = headerSize + sizeof (BitElem)
      * bitElemsAllocSize(bitsPerSeqpos * numMapEntries(
                            seqLen, mapIntervalLog2));
    mapName = str_clone(projectName);
    {
      char buf[1 + 4 + 3];
      snprintf(buf, sizeof (buf), ".%ucxm", (unsigned)mapIntervalLog2);
      str_append_cstr(mapName, buf);
      if (createMapFile)
      {
        /* write header information, for reference to verify noone
         * toyed with the file name */
        BitElem headerBuf[headerBitElems];
        if (!(mapFile = fa_fopen(str_get(mapName), "w+b")))
          break;
        bsStoreUInt16(headerBuf, 0, HEADER_ENTRY_BITS, mapIntervalLog2);
        bsStoreUInt16(headerBuf, HEADER_ENTRY_BITS, HEADER_ENTRY_BITS,
                      bitsPerSeqpos);
        if (fwrite(headerBuf,  sizeof (headerBuf), 1, mapFile) != 1)
          break;
        if (fseeko(mapFile, mapSize - 1, SEEK_SET))
          break;
        /* write one byte so mmap works for full file */
        if (fwrite(&buf, 1, 1, mapFile) != 1)
          break;
        if (fflush(mapFile) == EOF)
          break;
      }
      else
      {
        /* read header information, for reference to verify noone
         * toyed with the file name */
        BitElem headerBuf[headerBitElems];
        if (!(mapFile = fa_fopen(str_get(mapName), "rb")))
          break;
        if (fread(headerBuf,  sizeof (headerBuf), 1, mapFile) != 1)
          break;
        if (bsGetUInt16(headerBuf, 0, HEADER_ENTRY_BITS) != mapIntervalLog2
            || (bsGetUInt16(headerBuf, HEADER_ENTRY_BITS, HEADER_ENTRY_BITS)
                != bitsPerSeqpos))
        {
          fprintf(stderr, "error: context map file %s contains corrupted "
                  "data.\n", str_get(mapName));
          break;
        }
      }
    }
    mapMap = fa_mmap_generic_fd(fileno(mapFile), mapSize, 0, createMapFile,
                                false);
    newBWTSeqCR->revMap = (newBWTSeqCR->revMapMMap = mapMap) + headerSize;
  } while (0);
  if (mapName) str_delete(mapName);
  if (mapFile) fa_xfclose(mapFile);
  return mapMap != NULL;
}

extern BWTSeqContextRetriever *
BWTSeqCRLoad(const BWTSeq *bwtSeq, const Str *projectName,
             short mapIntervalLog2)
{
  Seqpos seqLen;
  unsigned short bitsPerSeqpos;
  BWTSeqContextRetriever *newBWTSeqCR;
  assert(bwtSeq && projectName);
  seqLen = BWTSeqLength(bwtSeq);
  bitsPerSeqpos = requiredSeqposBits(seqLen - 1);
  newBWTSeqCR = ma_malloc(sizeof (*newBWTSeqCR));
  newBWTSeqCR->bitsPerSeqpos = bitsPerSeqpos;
  newBWTSeqCR->bwtSeq = bwtSeq;
  if (mapIntervalLog2 != CTX_MAP_ILOG_AUTOSIZE)
  {
    if (!BWTSeqCRMapOpen(mapIntervalLog2, bitsPerSeqpos, seqLen,
                         projectName, false, newBWTSeqCR))
    {
      ma_free(newBWTSeqCR);
      return NULL;
    }
  }
  else
  {
    short mapIntervalLog2Max = bitsPerSeqpos;
    mapIntervalLog2 = 0;
    while (!BWTSeqCRMapOpen(mapIntervalLog2, bitsPerSeqpos, seqLen,
                            projectName, false, newBWTSeqCR)
           && mapIntervalLog2 < mapIntervalLog2Max)
      ++mapIntervalLog2;
    if (!newBWTSeqCR->revMap)
    {
      ma_free(newBWTSeqCR);
      return NULL;
    }
  }
  newBWTSeqCR->mapIntervalLog2 = mapIntervalLog2;
  newBWTSeqCR->mapInterval = 1 << mapIntervalLog2;
  newBWTSeqCR->mapMask = newBWTSeqCR->mapInterval - 1;
  return newBWTSeqCR;
}

extern void
deleteBWTSeqCR(BWTSeqContextRetriever *bwtSeqCR)
{
  assert(bwtSeqCR);
  fa_xmunmap(bwtSeqCR->revMapMMap);
  ma_free(bwtSeqCR);
}

extern void
BWTSeqCRAccessSubseq(const BWTSeqContextRetriever *bwtSeqCR,
                     Seqpos start, size_t len, Symbol subseq[])
{
  struct SeqMark currentPos;
  struct extBitsRetrieval extBits;
  const BWTSeq *bwtSeq;
  initExtBitsRetrieval(&extBits);
  assert(bwtSeqCR);
  assert(start < BWTSeqLength(bwtSeqCR->bwtSeq));
  bwtSeq = bwtSeqCR->bwtSeq;
  {
    Seqpos end = start + len - 1;
    currentPos = BWTSeqCRNextMark(bwtSeqCR, end);
    assert(currentPos.textPos >= end);
    while (currentPos.textPos > end)
    {
      currentPos.bwtPos = BWTSeqLFMap(bwtSeq, currentPos.bwtPos, &extBits);
      --currentPos.textPos;
    }
  }
  {
    Seqpos i = len;
    while (i)
    {
      --i;
      subseq[i] = BWTSeqGetSym(bwtSeq, currentPos.bwtPos);
      currentPos.bwtPos = BWTSeqLFMap(bwtSeq, currentPos.bwtPos, &extBits);
#ifndef NDEBUG
      --currentPos.textPos;
#endif
    }
  }
  assert(currentPos.textPos + 1 == start);
}
