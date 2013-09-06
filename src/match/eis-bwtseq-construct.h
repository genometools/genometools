/*
  Copyright (c) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>

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

#ifndef EIS_BWTSEQ_CONSTRUCT_H
#define EIS_BWTSEQ_CONSTRUCT_H

/**
 * \file eis-bwtseq-construct.h
 * Interface definitions for construction of an indexed representation of the
 * BWT of a sequence as presented by Manzini and Ferragina (Compressed
 * Representations of Sequences and Full-Text Indexes, 2006)
 * \author Thomas Jahns <Thomas.Jahns@gmx.net>
 */

#include "match/sarr-def.h"
#include "match/eis-bwtseq.h"
#include "match/eis-suffixerator-interface.h"
#include "match/eis-suffixarray-interface.h"

/**
 * \brief Loads (or creates if necessary) an encoded indexed sequence
 * object of the BWT transform.
 * @param params a struct holding parameter information for index construction
 * @param sa Suffixarray data structure to build BWT index from
 * @param totalLen length of sorted sequence (including terminator)
 * @param err genometools error object reference
 * @return reference to new BWT sequence object
 */
BWTSeq *
gt_availBWTSeqFromSA(const struct bwtParam *params, Suffixarray *sa,
                  GtUword totalLen, GtError *err);

/**
 * \brief Loads an encoded indexed sequence object of the
 * BWT transform.
 * @param params a struct holding parameter information for index construction
 * @param sa Suffixarray data structure to build BWT index from
 * @param totalLen length of BWT sequence (including terminator symbol)
 * @param err genometools error object reference
 * @return reference to new BWT sequence object
 */
BWTSeq *
gt_loadBWTSeqForSA(const char *projectName, enum seqBaseEncoding encType,
                int BWTOptFlags, const GtAlphabet *gtalphabet,
                GtError *err);

/**
 * \brief Creates an encoded indexed sequence object of the BWT transform.
 * @param params a struct holding parameter information for index construction
 * @param si Suffixerator interface to read data for BWT index from
 * @param totalLen length of BWT sequence (including terminator symbol)
 * @param err genometools reference for core functions
 * @return reference to new BWT sequence object
 */
BWTSeq *
gt_createBWTSeqFromSfxI(const struct bwtParam *params, sfxInterface *si,
                     GtError *err);

/**
 * \brief Creates or loads an encoded indexed sequence object of the
 * BWT transform.
 * @param params a struct holding parameter information for index construction
 * @param sa suffix array to read data for BWT index from
 * @param totalLen length of BWT sequence (including terminator symbol)
 * @param err genometools reference for core functions
 * @return reference to new BWT sequence object
 */
BWTSeq *
gt_createBWTSeqFromSA(const struct bwtParam *params, Suffixarray *sa,
                   GtUword totalLen, GtError *err);

/**
 * \brief Creates or loads an encoded indexed sequence object of the
 * BWT transform.
 * @param params a struct holding parameter information for index construction
 * @param sai interface to suffix array to read data for BWT index from
 * @param totalLen length of BWT sequence (including terminator symbol)
 * @param err genometools reference for core functions
 * @return reference to new BWT sequence object
 */
BWTSeq *
gt_createBWTSeqFromSAI(const struct bwtParam *params,
                    SuffixarrayFileInterface *sai,
                    GtError *err);

#endif
