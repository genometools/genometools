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
#include "core/queue.h"
#include "core/str_array.h"

/* the GtSequenceBuffer interface */
typedef struct GtSequenceBuffer GtSequenceBuffer;
typedef struct GtSequenceBufferClass GtSequenceBufferClass;

/* Increases the reference count of the <GtSequenceBuffer>. */
GtSequenceBuffer*  gt_sequence_buffer_ref(GtSequenceBuffer*);

/* Creates a new <GtSequenceBuffer>, choosing the appropriate type by looking
   at the first input file. All files must be of the same type.
   If NULL is returned, an error occurred. */
GtSequenceBuffer*  gt_sequence_buffer_new_guess_type(GtStrArray*, GtError*);

/* Fetches next character from <GtSequenceBuffer>. */
int           gt_sequence_buffer_next(GtSequenceBuffer*, unsigned char*,
                                      GtError*);

/* Advances the sequence window in the <GtSequenceBuffer> by OUTBUFSIZE. */
int           gt_sequence_buffer_advance(GtSequenceBuffer*, GtError*);

/* Returns the index of the currently read sequence file. */
unsigned long gt_sequence_buffer_get_file_index(GtSequenceBuffer*);

/* Assigns a symbol map to the sequence iterator to transform sequences with.
   Set to NULL to disable alphabet transformation (default). */
void          gt_sequence_buffer_set_symbolmap(GtSequenceBuffer*,
                                               const unsigned char*);

/* Assigns an array of Filelengthvalue structs to the sequence iterator which
   is filled during iteration. Note that the length of the array must equal the
   number of sequence files traversed.
   Set to NULL to disable Filelengthvalue counting (default). */
void          gt_sequence_buffer_set_filelengthtab(GtSequenceBuffer*,
                                                   Filelengthvalues*);

void          gt_sequence_buffer_set_desc_queue(GtSequenceBuffer *si,
                                                GtQueue *dq);

/* Assigns an array of sizeof (char) length which counts the occurrences of each
   alphabet character in the read sequence.
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

int           gt_sequence_buffer_unit_test(GtError*);

void          gt_sequence_buffer_delete(GtSequenceBuffer*);

#endif
