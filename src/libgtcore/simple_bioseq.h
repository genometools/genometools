/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#ifndef SIMPLE_BIOSEQ_H
#define SIMPLE_BIOSEQ_H

/*
   This simple class allows to parse fasta files and store them in a fashion
   which is conveniently accessible.

   It implements a subset of the behavior of the Bioseq class which can be found
   in GenomeTools, but with much fewer dependencies and considerably reduced
   funcionality. It serves mainly educational purposes and should be easily
   understandable.
*/
typedef struct SimpleBioseq SimpleBioseq;

SimpleBioseq* simple_bioseq_new(const char *fasta_file);
void          simple_bioseq_delete(SimpleBioseq*);
const char*   simple_bioseq_get_description(SimpleBioseq*, unsigned long);
const char*   simple_bioseq_get_sequence(SimpleBioseq*, unsigned long index);
unsigned long simple_bioseq_get_sequence_length(SimpleBioseq*, unsigned long);
unsigned long simple_bioseq_number_of_sequences(SimpleBioseq*);

#endif
