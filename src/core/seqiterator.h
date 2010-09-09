/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007-2010 Center for Bioinformatics, University of Hamburg

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

#ifndef SEQITERATOR_H
#define SEQITERATOR_H

#include <inttypes.h>
#include "core/queue.h"
#include "core/str_array.h"
#include "core/sequence_buffer.h"
#include "core/types_api.h"

typedef struct GtSeqIterator GtSeqIterator;
typedef struct GtSeqIteratorClass GtSeqIteratorClass;

/* Sets a symbol map for the <GtSeqIterator>.
   If a <symbolmap> is given, all read in sequences are transformed with it.
   Set to NULL to disable alphabet transformation. */
void           gt_seqiterator_set_symbolmap(GtSeqIterator*,
                                            const GtUchar *symbolmap);

/* If set to <true>, sequences and descriptions are processed (otherwise
   only the descriptions). By default, sequences are processed. */
void           gt_seqiterator_set_sequence_output(GtSeqIterator*, bool);

/* Get next <sequence> (of length <len>) and <description> from <seqit>.
   The caller is responsible to free the received <description>.
   Returns 1, if another sequence could be parsed. 0, if all given sequence
   files are exhausted. And -1, if an error occured (<err> is set
   accordingly). */
int            gt_seqiterator_next(GtSeqIterator *seqit,
                                   const GtUchar **sequence,
                                   unsigned long *len,
                                   char **description, GtError*);

/* Returns a pointer to the current total number of read characters. */
const unsigned
long long*     gt_seqiterator_getcurrentcounter(GtSeqIterator*,
                                                unsigned long long);

/* Returns TRUE if <seqit> supports setting of a quality buffer. */
bool           gt_seqiterator_has_qualities(GtSeqIterator *seqit);

/* Turns on reporting of sequence qualities to the location
   pointed to by <qualities>. That pointer will be set to a string containing
   the quality data (which must then be processed into scores). */
void           gt_seqiterator_set_quality_buffer(GtSeqIterator *seqit,
                                                 const GtUchar **qualities);

/* Deletes <seqit> and frees associated memory. */
void           gt_seqiterator_delete(GtSeqIterator *seqit);

#endif
