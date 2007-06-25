/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtext/reverse.h>
#include <libgtext/splicedseq.h>

struct Splicedseq {
  Str *splicedseq;
  Array *positionmapping;
  bool forward;
};

Splicedseq* splicedseq_new(Env *env)
{
  Splicedseq *ss = env_ma_malloc(env, sizeof (Splicedseq));
  ss->splicedseq = str_new(env);
  ss->positionmapping = array_new(sizeof (unsigned long), env);
  ss->forward = true;
  return ss;
}

void splicedseq_add(Splicedseq *ss, unsigned long start, unsigned long end,
                    const char *original_sequence, Env *env)
{
  unsigned long i;
  assert(ss && start <= end && original_sequence);
  str_append_cstr_nt(ss->splicedseq, original_sequence + start,
                     end - start + 1, env);
  /* make sure elemnts are added in ascending order */
  assert(!array_size(ss->positionmapping) ||
         start > *(unsigned long*) array_get_last(ss->positionmapping));
  for (i = start; i <= end; i++)
    array_add(ss->positionmapping, i, env);
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

int splicedseq_reverse(Splicedseq *ss, Env *env)
{
  int had_err;
  assert(ss);
  had_err = reverse_complement(str_get(ss->splicedseq),
                               str_length(ss->splicedseq), env);
  if (!had_err) {
    array_reverse(ss->positionmapping, env);
    ss->forward = !ss->forward;
  }
  return had_err;
}

void splicedseq_reset(Splicedseq *ss)
{
  assert(ss);
  str_reset(ss->splicedseq);
  array_reset(ss->positionmapping);
  ss->forward = true;
}

static int check_splicedseq(Splicedseq *ss, Env *env)
{                       /*0123456789*/
  static char *origseq = "aaccaagtga", *splicedseq = "ccgtg";
  int had_err = 0;
  env_error_check(env);
  splicedseq_add(ss, 2, 3, origseq, env);
  splicedseq_add(ss, 6, 8, origseq, env);
  ensure(had_err, strcmp(splicedseq_get(ss), splicedseq) == 0);
  ensure(had_err, !splicedseq_pos_is_border(ss, 0));
  ensure(had_err,  splicedseq_pos_is_border(ss, 1));
  ensure(had_err, !splicedseq_pos_is_border(ss, 2));
  ensure(had_err, !splicedseq_pos_is_border(ss, 3));
  ensure(had_err, !splicedseq_pos_is_border(ss, 4));
  return had_err;
}

int splicedseq_unit_test(Env *env)
{
  Splicedseq *ss;
  int had_err = 0;
  env_error_check(env);
  ss = splicedseq_new(env);
  had_err = check_splicedseq(ss, env);
  if (!had_err) {
    splicedseq_reset(ss);
    had_err = check_splicedseq(ss, env);
  }
  splicedseq_delete(ss, env);
  return had_err;
}

void splicedseq_delete(Splicedseq *ss, Env *env)
{
  if (!ss) return;
  str_delete(ss->splicedseq, env);
  array_delete(ss->positionmapping, env);
  env_ma_free(ss, env);
}
