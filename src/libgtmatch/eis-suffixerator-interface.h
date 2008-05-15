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
/**
 * \file eis-suffixerator-interface.h
 * \brief Methods to call suffixerator functions through one object,
 * but have the same data available to multiple listeners.
 * Conforms to the abstract interface defined in
 * eis-sa-common.h for suffix array class objects (type SASeqSrc).
 * \author Thomas Jahns <Thomas.Jahns@gmx.net>
 */
#ifndef EIS_SUFFIXERATOR_INTERFACE_H
#define EIS_SUFFIXERATOR_INTERFACE_H

#include "libgtcore/error.h"
#include "libgtmatch/sfx-suffixer.h"
#include "libgtmatch/eis-mrangealphabet.h"
#include "libgtmatch/eis-seqdatasrc.h"
#include "libgtmatch/eis-sa-common.h"

/**
 * opaque interface layer to retrieve arbitrary length portions of
 * the suffixarray
 */
typedef struct sfxInterface sfxInterface;

/**
 * @brief Create suffixerator interface object.
 *
 * @param so options for calling suffixerator
 * @param encseq object holding the sequences for suffixerator to
 * operate on
 * @param specialcharinfo
 * @param numofsequences number of sequences concatenated in encseq
 * @param mtime timing device
 * @param length length of concatenated sequences plus terminator and
 * separators
 * @param alpha alphabet to use
 * @param characterdistribution counts for all characters, used to
 * generate statistics
 * @param verbosity used as argument of showverbose
 * @param err
 * @return interface object reference
 */
extern sfxInterface *
newSfxInterface(Readmode readmode,
                unsigned int prefixlength,
                unsigned int numofparts,
                const Sfxstrategy *sfxstrategy,
                const Encodedsequence *encseq,
                const Specialcharinfo *specialcharinfo,
                unsigned long numofsequences,
                Measuretime *mtime,
                Seqpos length,
                const Alphabet *alpha,
                const unsigned long *characterdistribution,
                Verboseinfo *verbosity,
                Error *err);

/**
 * @brief Create suffixerator interface object with requestors already
 * in place.
 *
 * @param so options for calling suffixerator
 * @param numReaders number of readers to register at creation
 * @param requests reference of array with request type for each
 * reader
 * @param ids array to which, for each reader, one id is written (same
 * sequence as corresponding requests)
 * @param encseq object holding the sequences for suffixerator to
 * operate on
 * @param specialcharinfo
 * @param numofsequences number of sequences concatenated in encseq
 * @param mtime timing device
 * @param length length of concatenated sequences plus terminator and
 * separators
 * @param alpha alphabet to use
 * @param characterdistribution counts for all characters, used to
 * generate statistics
 * @param verbosity used as argument of showverbose
 * @param err
 * @return interface object reference
 */
extern sfxInterface *
newSfxInterfaceWithReaders(Readmode readmode,
                           unsigned int prefixlength,
                           unsigned int numofparts,
                           const Sfxstrategy *sfxstrategy,
                           size_t numReaders,
                           enum sfxDataRequest readerRequests[],
                           SeqDataReader readers[],
                           const Encodedsequence *encseq,
                           const Specialcharinfo *specialcharinfo,
                           unsigned long numofsequences,
                           Measuretime *mtime,
                           Seqpos length,
                           const Alphabet *alpha,
                           const unsigned long *characterdistribution,
                           Verboseinfo *verbosity,
                           Error *err);

/**
 * @brief get Sfxiterator from SfxInterface
 *
 * @param iface to take Sfxiterator from
 *
 */
extern const Sfxiterator *
SfxInterface2Sfxiterator(const sfxInterface *iface);

/**
 * @brief Deallocate resources of suffixerator interface object.
 *
 * @param iface object to delete
 */
extern void
deleteSfxInterface(sfxInterface *iface);

/**
 * @brief Dynamically cast to super class.
 *
 * @param sfxi reference to suffixerator interface object
 * @return reference of object of base class
 */
extern SASeqSrc *
SfxI2SASS(sfxInterface *sfxi);

/**
 * \brief Constructs multiple range alphabet for sequence sorted by
 * suffixerator (i.e. alphabet includes separator symbol).
 *
 * @param si reference of interface to suffixerator
 * @return reference of newly created alphabet object
 */
extern MRAEnc *
SfxINewMRAEnc(const sfxInterface *si);

/**
 * \brief Get reference for alphabet used to encode original sequence
 * object.
 *
 * @param si reference of interface to suffixerator
 * @return reference of alphabet object
 */
extern const Alphabet *
SfxIGetAlphabet(const sfxInterface *si);

/**
 * \brief Get reference for original sequence object.
 *
 * @param si reference of interface to suffixerator
 * @return reference of sequence object
 */
extern const Encodedsequence *
SfxIGetEncSeq(const sfxInterface *si);

/**
 * @brief Get read mode used for suffix sorting.
 * @param si suffixerator interface object reference
 * @return read mode
 */
extern Readmode
SfxIGetReadmode(const sfxInterface *si);

/**
 * \brief Get original sequence substring.
 *
 * @param si reference of interface to suffixerator
 * @param dest store read symbols here
 * @param pos position to start reading at
 * @param len number of symbols to read
 * @return number of symbols actually read
 */
extern size_t
SfxIGetOrigSeq(const void *si, Symbol *dest, Seqpos pos, size_t len);

/**
 * \brief Query original sequence for statistics.
 *
 * @param si reference of interface to suffixerator
 * @return reference of struct holding statistics (symbol counts)
 */
extern const struct seqStats *
SfxIGetSeqStats(const sfxInterface *si);

/**
 * @brief Query length @f$l@f$ of sequence sorted by suffixerator, including
 * the terminator and separator symbols (i.e.
 * @f[ l = \sum_{i=1}^n \left(|s_i| + 1\right)@f]
 * ).
 *
 * @param si reference of interface to suffixerator
 * @return length of sequence
 */
extern Seqpos
SfxIGetLength(const sfxInterface *si);

/**
 * @brief Query position of suffix starting at position 0, can be
 * undefined if not yet encountered.
 *
 * @param si reference of interface to suffixerator
 * @return tuple of boolean (position is known) and position (if
 * known) or undefined value.
 */
extern DefinedSeqpos
SfxIGetRot0Pos(const struct sfxInterface *si);

/**
 * @return >0 on success, 0 on error
 */
extern SeqDataReader
SfxIRegisterReader(sfxInterface *iface, enum sfxDataRequest request);

#if 0
/**
 * \brief Reads portion of the BWT string produced by suffixerator.
 *
 * @param iface
 * @param id value returned by corresponding SfxIRegisterReader call
 * @param len number of symbols to read
 * @param dest store read symbols here
 * @return number of symbols read (less than len implies end of file)
 */
extern size_t
readSfxIBWTRange(sfxInterface *iface, listenerID id, size_t len, Uchar *dest);

/**
 * @return actual number of symbols read
 */
extern size_t
readSfxILCPRange(sfxInterface *iface, listenerID id, size_t len, Seqpos *dest,
                 Error *err);

/**
 * @return actual number of symbols read
 */
extern size_t
readSfxISufTabRange(sfxInterface *iface, listenerID id, size_t len,
                    Seqpos *dest);
#endif

extern SeqDataReader
SfxIRegisterReader(sfxInterface *iface, enum sfxDataRequest request);

#endif
