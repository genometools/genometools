/*
  Copyright (c) 2013 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

#ifndef HPLSTORE_H
#define HPLSTORE_H

#include <stdint.h>
#include "core/encseq_api.h"
#include "core/file_api.h"

typedef struct GtHplstore GtHplstore;

/* The <GtHplstore> class provides a data structure for the storage
   of an homopolymer lengths array over a collection of sequences
   (such as sequencing reads). */

/* <GT_HPLSTORE_MAX> is the maximal homopolymer length allowed */
#define    GT_HPLSTORE_MAX (uint8_t)254U

/* <GT_HPLSTORE_RANGE> is the number of homopolymer values which
   may be represented*/
#define    GT_HPLSTORE_RANGE (uint8_t)(GT_HPLSTORE_MAX + 1U)

/* <GT_HPLSTORE_UNDEF> is a reserved value which is never stored in an
   <GtHplstore> and thus can be safely used to represent undefined values */
#define    GT_HPLSTORE_UNDEF (uint8_t)(GT_HPLSTORE_MAX + 1U)

/* create a new <GtHplstore> for a collection of sequences with
   <nofelements> symbols including the separators */
GtHplstore *gt_hplstore_new(GtUword nofelements);

/* deletes <hplstore> */
void       gt_hplstore_delete(GtHplstore *hplstore);

/* finalize <hplstore> (which from this point becomes read-only);
   the finalized hplstore contains <nofelements> values, which may be equal or
   less than the initial nofelements used for <gt_hplstore_new>;
   a <GtHplstore> may only be finalized once */
void       gt_hplstore_finalize(GtHplstore *hplstore,
                                GtUword nofelements);

/* sets the value at position <pos> of <hplstore> to <value>;
   this method may only be called after <gt_hplstore_finalize> */
void       gt_hplstore_set(GtHplstore *hplstore, GtUword pos,
                           uint8_t value);

/* gets the value at position <pos> of <hplstore>;
   this method may only be called before <gt_hplstore_finalize> */
uint8_t    gt_hplstore_get(GtHplstore *hplstore, GtUword pos);

/* copy <nofelements> values stored in the <GtHplstore> from position
   <from>in the array pointed by <hplengths>;
   the caller is responsible to allocate enough space in <hplengths>;
   this method may only be called after <gt_hplstore_finalize> */
void       gt_hplstore_get_range(const GtHplstore *hplstore, uint8_t *hplengths,
                                 GtUword from, GtUword nofelements);

/* show to <outfile> the sequence of <nofelements> homopolymers starting from
   homopolymer position <from>, obtained by combining the homopolymer
   lengths stored in <hplstore> with the symbols stored in <encseq>;
   this method may only be called after <gt_hplstore_finalize> */
void       gt_hplstore_show_decoded_sequence(GtFile *outfile,
                                             const GtHplstore *hplstore,
                                             const GtEncseq *encseq,
                                             GtUword from,
                                             GtUword nofelements);

/* more efficient variant of <gt_hplstore_show_decoded_sequence> in which the
   length component is provided as an array such that obtained by the method
   <gt_hplstore_get_range>; <encseq_from> is the offset in the
   encseq; the homopolymer lengths array is of length <nofelements> */
void       gt_hplstore_show_decoded_sequence_using_hplengths(
                                             GtFile *outfile,
                                             const uint8_t *hplengths,
                                             const GtEncseq *encseq,
                                             GtUword encseq_from,
                                             GtUword nofelements);

/* save a finalized <GtHplstore> to the file <out_fp>;
   the number of values must be stored separately elsewhere;
   this method may only be called after <gt_hplstore_finalize> */
void       gt_hplstore_save(const GtHplstore *hplstore, FILE *out_fp);

/* load a finalized <GtHplstore> containing <nofelements> values
   from the file <in_fp> */
GtHplstore *gt_hplstore_load(FILE *in_fp, GtUword nofelements);

#endif
