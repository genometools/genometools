/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "array.h"
#include "ensure.h"
#include "reverse.h"
#include "splicedseq.h"
#include "str.h"
#include "xansi.h"

struct Splicedseq {
  Str *splicedseq;
  Array *positionmapping;
  bool forward;
};

Splicedseq* splicedseq_new(void)
{
  Splicedseq *ss = xmalloc(sizeof(Splicedseq));
  ss->splicedseq = str_new();
  ss->positionmapping = array_new(sizeof(unsigned long));
  ss->forward = true;
  return ss;
}

void splicedseq_add(Splicedseq *ss, unsigned long start, unsigned long end,
                    const char *original_sequence)
{
  unsigned long i;
  assert(ss && start <= end && original_sequence);
  str_append_cstr_nt(ss->splicedseq, original_sequence + start,
                     end - start + 1);
  /* make sure elemnts are added in ascending order */
  assert(!array_size(ss->positionmapping) ||
         start > *(unsigned long*) array_get_last(ss->positionmapping));
  for (i = start; i <= end; i++)
    array_add(ss->positionmapping, i);
}

char* splicedseq_get(const Splicedseq *ss)
{
  return str_get(ss->splicedseq);
}

bool splicedseq_pos_is_border(const Splicedseq *ss, unsigned long pos)
{
  assert(ss && str_length(ss->splicedseq) == array_size(ss->positionmapping));
  assert(pos < str_length(ss->splicedseq)); /* legal position */
  if (ss->forward && pos + 1 < array_size(ss->positionmapping) &&
      *(unsigned long*) array_get(ss->positionmapping, pos) + 1 !=
      *(unsigned long*) array_get(ss->positionmapping, pos+1)) {
    return true;
  }
  if (!ss->forward && pos &&
      *(unsigned long*) array_get(ss->positionmapping, pos-1) - 1 !=
      *(unsigned long*) array_get(ss->positionmapping, pos)) {
    return true;
  }
  return false;
}

unsigned long splicedseq_map(const Splicedseq *ss, unsigned long pos)
{
  assert(ss && str_length(ss->splicedseq) == array_size(ss->positionmapping));
  assert(pos < str_length(ss->splicedseq)); /* legal position */
  return *(unsigned long*) array_get(ss->positionmapping, pos);
}

unsigned long splicedseq_length(const Splicedseq *ss)
{
  assert(ss);
  return str_length(ss->splicedseq);
}

void splicedseq_reverse(Splicedseq *ss)
{
  /*unsigned long i;*/
  assert(ss);
  reverse_complement(str_get(ss->splicedseq), str_length(ss->splicedseq));
#if 0
  for (i = 0; i < array_size(ss->positionmapping); i++) {
    printf("%lu ", *(unsigned long*) array_get(ss->positionmapping, i));
  }
  putchar('\n');
#endif
  array_reverse(ss->positionmapping);
#if 0
  for (i = 0; i < array_size(ss->positionmapping); i++) {
    printf("%lu ", *(unsigned long*) array_get(ss->positionmapping, i));
  }
  putchar('\n');
#endif
  ss->forward = !ss->forward;
}

void splicedseq_reset(Splicedseq *ss)
{
  assert(ss);
  str_reset(ss->splicedseq);
  array_set_size(ss->positionmapping, 0);
  ss->forward = true;
}

static void check_splicedseq(Splicedseq *ss)
{                       /*0123456789*/
  static char *origseq = "aaccaagtga", *splicedseq = "ccgtg";

  splicedseq_add(ss, 2, 3, origseq);
  splicedseq_add(ss, 6, 8, origseq);
  ensure(strcmp(splicedseq_get(ss), splicedseq) == 0);
  ensure(!splicedseq_pos_is_border(ss, 0));
  ensure( splicedseq_pos_is_border(ss, 1));
  ensure(!splicedseq_pos_is_border(ss, 2));
  ensure(!splicedseq_pos_is_border(ss, 3));
  ensure(!splicedseq_pos_is_border(ss, 4));
}

int splicedseq_unit_test(void)
{
  Splicedseq *ss = splicedseq_new();
  check_splicedseq(ss);
  splicedseq_reset(ss);
  check_splicedseq(ss);
  splicedseq_free(ss);
  return EXIT_SUCCESS;
}

void splicedseq_free(Splicedseq *ss)
{
  if (!ss) return;
  str_free(ss->splicedseq);
  array_free(ss->positionmapping);
  free(ss);
}
