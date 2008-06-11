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

#include "libgtcore/cstr.h"
#include "libgtcore/ma.h"
#include "libgtexercise/sspliced_alignment.h"

struct SSplicedAlignment {
  char *id;
  bool forward;
  Array *exons; /* the exon ranges */
};

SSplicedAlignment* sspliced_alignment_new(const char *id, bool forward)
{
  SSplicedAlignment *sa;
  assert(id);
  sa = ma_malloc(sizeof *sa);
  sa->id = cstr_dup(id);
  sa->forward = forward;
  sa->exons = array_new(sizeof (Range));
  return sa;
}

void sspliced_alignment_delete(SSplicedAlignment *sa)
{
  if (!sa) return;
  array_delete(sa->exons);
  ma_free(sa->id);
  ma_free(sa);
}

bool sspliced_alignment_is_forward(const SSplicedAlignment *sa)
{
  assert(sa);
  return sa->forward;
}

void sspliced_alignment_add_exon(SSplicedAlignment *sa, Range exon)
{
  assert(sa);
  array_add(sa->exons, exon);
}

unsigned long sspliced_alignment_num_of_exons(const SSplicedAlignment *sa)
{
  assert(sa);
  return array_size(sa->exons);
}

Range sspliced_alignment_get_exon(const SSplicedAlignment *sa,
                                  unsigned long exon_number)
{
  assert(sa);
  return *(Range*) array_get(sa->exons, exon_number);
}

Range sspliced_alignment_genomic_range(const SSplicedAlignment *sa)
{
  Range range;
  assert(sa);
  assert(array_size(sa->exons));
  range.start = ((Range*) array_get_first(sa->exons))->start;
  range.end   = ((Range*) array_get_last(sa->exons))->end;
  return range;
}
