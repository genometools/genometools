/*
  Copyright (c) 2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#ifndef SEQUENCE_BUFFER_H
#define SEQUENCE_BUFFER_H

#include "core/chardef.h"
#include "core/error_api.h"
#include "core/filelengthvalues.h"
#include "core/desc_buffer.h"
#include "core/str_array.h"

/* A <GtSequenceBuffer> represents a group of sequence files of a certain type.
   These files are parsed on-the-fly and the sequences and descriptions
   contained in them are made available.
   Sequences can be read character-wise using gt_sequence_buffer_next(), with
   SEPARATOR symbols in between (see chardef.h).
   Note that the <GtSequenceBuffer> is a rather low-level tool for efficient
   sequence access. For simple access to whole sequences, use the
   <GtSeqIterator> class. */
typedef struct GtSequenceBuffer GtSequenceBuffer;
typedef struct GtSequenceBufferClass GtSequenceBufferClass;

/* Increases the reference count of the <GtSequenceBuffer>. */
GtSequenceBuffer*  gt_sequence_buffer_ref(GtSequenceBuffer*);

/* Creates a new <GtSequenceBuffer>, choosing the appropriate type by looking
   at the first input file. All files must be of the same type.
   If NULL is returned, an error occurred. */
GtSequenceBuffer*  gt_sequence_buffer_new_guess_type(const GtStrArray*,
                                                     GtError*);

/* Fetches next character from <GtSequenceBuffer>.
   Returns 1 if a new character could be read, 0 if all files are exhausted, or
   -1 on error (see the <GtError> object for details). */
int           gt_sequence_buffer_next(GtSequenceBuffer*, GtUchar*, GtError*);

/* Fetches next character from <GtSequenceBuffer>.
   This method also always delivers the original character at the current
   reading position, regardless of symbol mappings that may apply.
   Returns 1 if a new character could be read, 0 if all files are exhausted, or
   -1 on error (see the <GtError> object for details). */
int           gt_sequence_buffer_next_with_original(GtSequenceBuffer*,
                                                    GtUchar *val, char *orig,
                                                    GtError*);

/* Returns the index of the currently read sequence file in the input file
   <GtStrArray>. */
unsigned long gt_sequence_buffer_get_file_index(GtSequenceBuffer*);

/* Assigns a symbol map to the sequence iterator to transform sequences with.
   Set to NULL to disable alphabet transformation (default). */
void          gt_sequence_buffer_set_symbolmap(GtSequenceBuffer*,
                                               const GtUchar *);

/* Assigns an array of Filelengthvalue structs to the sequence iterator. This
   is filled during iteration. Note that the length of the array must equal the
   number of sequence files traversed.
   Set to NULL to disable Filelengthvalue counting (default). */
void          gt_sequence_buffer_set_filelengthtab(GtSequenceBuffer*,
                                                   GtFilelengthvalues*);

/* Assigns a <GtDescBuffer> in which for each sequence file, the respective
   description string is written. If this is not set, or set to NULL,
   then descriptions are ignored in the input files. */
void          gt_sequence_buffer_set_desc_buffer(GtSequenceBuffer *si,
                                                 GtDescBuffer *db);

/* Assigns an array which counts the occurrences of each alphabet character in
   the read sequence. It must have at least as many elements as the number of
   characters in the expected alphabet.
   Set to NULL to disable character distribution counting (default). */
void          gt_sequence_buffer_set_chardisttab(GtSequenceBuffer*,
                                                 unsigned long*);

/* Returns the length of the last processed continuous stretch of special
   characters (wildcards or separators, see chardef.h). */
uint64_t      gt_sequence_buffer_get_lastspeciallength(const GtSequenceBuffer*);

/* Returns a pointer to a memory location holding the number of characters
   read altogether in this sequence of files. */
const unsigned long long*
              gt_sequence_buffer_get_counter(const GtSequenceBuffer *si);

void          gt_sequence_buffer_delete(GtSequenceBuffer*);

int           gt_sequence_buffer_unit_test(GtError*);

#endif
