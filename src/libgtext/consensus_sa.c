/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg

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

#include <assert.h>
#include "libgtcore/bittab.h"
#include "libgtcore/undef.h"
#include "libgtext/consensus_sa.h"

typedef struct {
  const void *set_of_sas;
  unsigned long number_of_sas;
  size_t size_of_sa;
  GetGenomicRangeFunc get_genomic_range;
  GetStrandFunc get_strand;
  GetExonsFunc get_exons;
  ProcessSpliceFormFunc process_splice_form;
  void *userdata;
} ConsensusSA;

#ifndef NDEBUG
static bool set_of_sas_is_sorted(const void *set_of_sas,
                                 unsigned long number_of_sas,
                                 size_t size_of_sa, GetGenomicRangeFunc
                                 get_genomic_range, Env *env)
{
  Range range_a, range_b;
  unsigned long max_end;
  const void *sa;

  assert(set_of_sas && number_of_sas && size_of_sa && get_genomic_range);

  /* get first range */
  range_a = get_genomic_range(set_of_sas);
  env_log_log(env, "-from %lu", range_a.start);
  max_end = range_a.end;

  /* loop over the other ranges */
  for (sa = set_of_sas + size_of_sa;
       sa < set_of_sas + number_of_sas * size_of_sa;
       sa += size_of_sa) {
    /* get second range */
    range_b = get_genomic_range(sa);

    /* compare */
    if (!(range_a.start <= range_b.start)) return false;

    /* make the second range the first */
    range_a = range_b;
    if (range_a.end > max_end)
      max_end = range_a.end;
  }
  env_log_log(env, "-to %lu", max_end);
  return true;
}
#endif

static Range extract_genomic_range(const ConsensusSA *csa, unsigned long sa)
{
  assert(csa && csa->set_of_sas && sa < csa->number_of_sas);
  return csa->get_genomic_range(csa->set_of_sas + csa->size_of_sa * sa);
}

static Strand extract_strand(const ConsensusSA *csa, unsigned long sa)
{
  Strand strand;
  assert(csa && csa->set_of_sas && sa < csa->number_of_sas);
  strand = csa->get_strand(csa->set_of_sas + csa->size_of_sa * sa);
  assert(strand == STRAND_FORWARD || strand == STRAND_REVERSE); /* XXX */
  return strand;
}

static void extract_exons(const ConsensusSA *csa, Array *exon_ranges,
                          unsigned long sa, Env *env)
{
  assert(csa && exon_ranges && csa->set_of_sas && sa < csa->number_of_sas);
  csa->get_exons(exon_ranges, csa->set_of_sas + csa->size_of_sa * sa, env);
  assert(array_size(exon_ranges));
  assert(ranges_are_sorted_and_do_not_overlap(exon_ranges));
}

static bool has_donor_site(Array *gene, unsigned long exon)
{
  assert(exon < array_size(gene));
  if (exon == array_size(gene) - 1)
    return false;
  return true;
}

static bool has_acceptor_site(Array *gene, unsigned long exon)
{
  assert(exon < array_size(gene));
  if (exon == 0)
    return false;
  return true;
}

static bool compatible(const ConsensusSA *csa,
                       unsigned long sa_1, unsigned long sa_2, Env *env)
{
  Array *exons_sa_1, *exons_sa_2;
  Range range_sa_1, range_sa_2;
  unsigned long i, j, num_of_exons_1, num_of_exons_2,
                start_1 = UNDEF_ULONG, start_2 = UNDEF_ULONG;
  bool start_values_set = false;
  const unsigned long fuzzlength = 0; /* XXX */

  assert(csa);

  /* check strands */
  if (extract_strand(csa, sa_1) != extract_strand(csa, sa_2)) return false;

  /* init */
  exons_sa_1 = array_new(sizeof (Range), env);
  exons_sa_2 = array_new(sizeof (Range), env);

  /* get ranges */
  range_sa_1 = extract_genomic_range(csa, sa_1);
  range_sa_2 = extract_genomic_range(csa, sa_2);

  if (!range_overlap(range_sa_1, range_sa_2)) {
    array_delete(exons_sa_1, env);
    array_delete(exons_sa_2, env);
    return false;
  }

  /* get exons */
  extract_exons(csa, exons_sa_1, sa_1, env);
  extract_exons(csa, exons_sa_2, sa_2, env);

  /* determine the first overlapping exon pair */
  i = 0;
  j = 0;
  num_of_exons_1 = array_size(exons_sa_1);
  num_of_exons_2 = array_size(exons_sa_2);
  while (i < num_of_exons_1 && j < num_of_exons_2) {
    if (range_overlap(*(Range*) array_get(exons_sa_1, i),
                      *(Range*) array_get(exons_sa_2, j))) {
      start_1 = i;
      start_2 = j;
      start_values_set = true;
      break;
    }
    if (((Range*) array_get(exons_sa_1, i))->start <
        ((Range*) array_get(exons_sa_2, j))->start) {
      i++;
    }
    else
      j++;
  }
  if (!start_values_set) {
    array_delete(exons_sa_1, env);
    array_delete(exons_sa_2, env);
    return false;
  }
  /* from now on the start values are set */

  assert(start_1 != UNDEF_ULONG && start_2 != UNDEF_ULONG);
  if (!(start_1 == 0 || start_2 == 0)) {
    /* no first segment could be maped */
    array_delete(exons_sa_1, env);
    array_delete(exons_sa_2, env);
    return false;
  }

  while (start_1 < num_of_exons_1 && start_2 < num_of_exons_2) {
    range_sa_1 = *((Range*) array_get(exons_sa_1, start_1));
    range_sa_2 = *((Range*) array_get(exons_sa_2, start_2));

    if (range_overlap(range_sa_1, range_sa_2)) {
      /* analyze acceptor sites */

      /* see if at least one exon has a acceptor site (on the left).
         in this case additional checks have to be performed.
         Otherwise, this exons are compatible (on the left side)
         because they overlap */
      if (has_acceptor_site(exons_sa_1, start_1) ||
          has_acceptor_site(exons_sa_2, start_2)) {
        if (has_acceptor_site(exons_sa_1, start_1) &&
            has_acceptor_site(exons_sa_2, start_2) &&
            range_sa_1.start!= range_sa_2.start) {
          /* the acceptor sites are different */
          array_delete(exons_sa_1, env);
          array_delete(exons_sa_2, env);
          return false;
        }
        else if (has_acceptor_site(exons_sa_1, start_1) &&
                 range_sa_2.start + fuzzlength < range_sa_1.start) {
          /* not within fuzzlength */
          array_delete(exons_sa_1, env);
          array_delete(exons_sa_2, env);
          return false;
        }
        else if (has_acceptor_site(exons_sa_2, start_2) &&
                 range_sa_1.start + fuzzlength < range_sa_2.start) {
          /* not within fuzzlength */
          array_delete(exons_sa_1, env);
          array_delete(exons_sa_2, env);
          return false;
        }
      }
      /* analyze donor sites */

      /* see if at least one exon has a donor site (on the right).
         in this case additional checks have to be performed.
         Otherwise, this exons are compatible (on the right side)
         because they overlap */
      if (has_donor_site(exons_sa_1, start_1) ||
          has_donor_site(exons_sa_2, start_2)) {
        if (has_donor_site(exons_sa_1, start_1) &&
            has_donor_site(exons_sa_2, start_2) &&
            range_sa_1.end != range_sa_2.end) {
          /* the donor sites are different */
          array_delete(exons_sa_1, env);
          array_delete(exons_sa_2, env);
          return false;
        }
        else if (has_donor_site(exons_sa_1, start_1) &&
                 range_sa_2.end - fuzzlength > range_sa_1.end) {
          /* not within fuzzlength */
          array_delete(exons_sa_1, env);
          array_delete(exons_sa_2, env);
          return false;
        }
        else if (has_donor_site(exons_sa_2, start_2) &&
                 range_sa_1.end - fuzzlength > range_sa_2.end) {
          /* not within fuzzlength */
          array_delete(exons_sa_1, env);
          array_delete(exons_sa_2, env);
          return false;
        }
      }
    }
    else {
      /* no overlap: two ordered segments do not overlap each other */
      array_delete(exons_sa_1, env);
      array_delete(exons_sa_2, env);
      return false;
    }
    start_1++;
    start_2++;
  }

  /* passed all tests */
  array_delete(exons_sa_1, env);
  array_delete(exons_sa_2, env);
  return true;
}

static bool contains(const ConsensusSA *csa,
                     unsigned long sa_1, unsigned long sa_2, Env *env)
{
  Range range_sa_1, range_sa_2;
  assert(csa);

  /* get ranges */
  range_sa_1 = extract_genomic_range(csa, sa_1);
  range_sa_2 = extract_genomic_range(csa, sa_2);

  if (range_contains(range_sa_1, range_sa_2) &&
      compatible(csa, sa_1, sa_2, env)) {
    return true;
  }
  return false;
}

static void compute_C(Bittab **C, const ConsensusSA *csa, Env *env)
{
  unsigned long sa, sa_1;
  assert(csa);
  for (sa = 0; sa < csa->number_of_sas; sa++) {
    for (sa_1 = 0; sa_1 < csa->number_of_sas; sa_1++) {
      if (contains(csa, sa, sa_1, env))
        bittab_set_bit(C[sa], sa_1);
    }
    assert(bittab_bit_is_set(C[sa], sa));
  }
}

static void compute_left_or_right(Bittab **left_or_right,
                                  const ConsensusSA *csa,
                                  bool (*cmp_func) (const ConsensusSA *csa,
                                                    unsigned long sa_1,
                                                    unsigned long sa_2),
                                                    Env *env)
{
  unsigned long sa, sa_1;
  assert(csa && left_or_right && *left_or_right);
  for (sa = 0; sa < csa->number_of_sas; sa++) {
    for (sa_1 = 0; sa_1 < csa->number_of_sas; sa_1++) {
      if (cmp_func(csa, sa, sa_1) && compatible(csa, sa, sa_1, env))
        bittab_set_bit(left_or_right[sa], sa_1);
    }
  }
}

static bool is_right_of(const ConsensusSA *csa,
                        unsigned long sa_1, unsigned long sa_2)
{
  Range range_sa_1, range_sa_2;
  assert(csa);
  range_sa_1 = extract_genomic_range(csa, sa_1);
  range_sa_2 = extract_genomic_range(csa, sa_2);
  if (range_sa_1.start  > range_sa_2.start && range_sa_1.end > range_sa_2.end)
    return true;
  return false;
}

static bool is_left_of(const ConsensusSA *csa,
                       unsigned long sa_1, unsigned long sa_2)
{
  Range range_sa_1, range_sa_2;
  assert(csa);
  range_sa_1 = extract_genomic_range(csa, sa_1);
  range_sa_2 = extract_genomic_range(csa, sa_2);
  if (range_sa_1.start  < range_sa_2.start && range_sa_1.end < range_sa_2.end)
    return true;
  return false;
}

static void compute_left(Bittab **left, const ConsensusSA *csa, Env *env)
{
  assert(csa);
  compute_left_or_right(left, csa, is_right_of, env);
}

static void compute_right(Bittab **right, const ConsensusSA *csa, Env *env)
{
  assert(csa);
  compute_left_or_right(right, csa, is_left_of, env);
}

static void compute_L(Bittab **L, Bittab **C, Bittab **left,
                      unsigned long number_of_sas, Env *env)
{
  unsigned long sa, sa_1, sa_2, sa_1_size = 0, sa_2_size;
  Bittab *tmpset = bittab_new(number_of_sas, env);

  for (sa = 0; sa < number_of_sas; sa++) {
    sa_1 = UNDEF_ULONG;

    if (!bittab_is_true(left[sa])) {
      /* bittab is empty */
      bittab_equal(L[sa], C[sa]);
    }
    else {
      for (sa_2  = bittab_get_first_bitnum(left[sa]);
           sa_2 != bittab_get_last_bitnum(left[sa]);
           sa_2  = bittab_get_next_bitnum(left[sa], sa_2)) {
        if (sa_1 == UNDEF_ULONG) {
          sa_1 = sa_2;

          bittab_or(tmpset, L[sa_1], C[sa]);
          sa_1_size = bittab_count_set_bits(tmpset);
        }
        else {
          bittab_or(tmpset, L[sa_2], C[sa]);
          sa_2_size = bittab_count_set_bits(tmpset);

          if (sa_2_size > sa_1_size) {
            sa_1      = sa_2;
            sa_1_size = sa_2_size;
          }
        }
      }

      assert(sa_1 != UNDEF_ULONG);
      bittab_or(L[sa], L[sa_1], C[sa]);
    }
  }
  bittab_delete(tmpset, env);
}

static void compute_R(Bittab **R, Bittab **C, Bittab **right,
                      unsigned long number_of_sas, Env *env)
{
  unsigned long sa_1, sa_2, sa_1_size = 0, sa_2_size;
  long sa;
  Bittab *tmpset = bittab_new(number_of_sas, env);

  for (sa = number_of_sas-1; sa >= 0; sa--) {
    sa_1 = UNDEF_ULONG;

    if (!bittab_is_true(right[sa])) {
      /* bittab is empty */
      bittab_equal(R[sa], C[sa]);
    }
    else {
      for (sa_2  = bittab_get_first_bitnum(right[sa]);
           sa_2 != bittab_get_last_bitnum(right[sa]);
           sa_2  = bittab_get_next_bitnum(right[sa], sa_2)) {
        if (sa_1 == UNDEF_ULONG) {
          sa_1 = sa_2;

          bittab_or(tmpset, R[sa_1], C[sa]);
          sa_1_size = bittab_count_set_bits(tmpset);
        }
        else {
          bittab_or(tmpset, R[sa_2], C[sa]);
          sa_2_size = bittab_count_set_bits(tmpset);

          if (sa_2_size > sa_1_size) {
            sa_1      = sa_2;
            sa_1_size = sa_2_size;
          }
        }
      }
      assert(sa_1 != UNDEF_ULONG);
      bittab_or(R[sa], R[sa_1], C[sa]);
    }
  }
  bittab_delete(tmpset, env);
}

#ifndef NDEBUG
static bool splice_form_is_valid(Bittab *SA_p, const ConsensusSA *csa,
                                 Env *env)
{
  Bittab *SA_p_complement; /* SA \ SA_p */
  unsigned long sa, sa_prime;
  bool incompatible_found, valid = true;

  SA_p_complement = bittab_new(csa->number_of_sas, env);
  bittab_complement(SA_p_complement, SA_p);

  for (sa_prime  = bittab_get_first_bitnum(SA_p_complement);
       sa_prime != bittab_get_last_bitnum(SA_p_complement);
       sa_prime  = bittab_get_next_bitnum(SA_p_complement, sa_prime)) {
    incompatible_found = false;
    for (sa  = bittab_get_first_bitnum(SA_p);
         sa != bittab_get_last_bitnum(SA_p);
         sa  = bittab_get_next_bitnum(SA_p, sa)) {
      if (!compatible(csa, sa, sa_prime, env)) {
        incompatible_found = true;
        break;
      }
    }
    if (!incompatible_found) { valid = false; break; }
  }
  bittab_delete(SA_p_complement, env);
  return valid;
}
#endif

static void compute_csas(ConsensusSA *csa, Env *env)
{
  unsigned long i, sa_i, sa_i_size = 0, sa_prime, sa_prime_size;
  Array *splice_form;
  Bittab **C, **left, **right, **L, **R, *U_i, *SA_i, *SA_prime;
#ifndef NDEBUG
  unsigned long u_i_size, u_i_minus_1_size;
  assert(csa && csa->set_of_sas);
#endif

  /* init sets */
  C     = env_ma_malloc(env, sizeof (Bittab*) * csa->number_of_sas);
  left  = env_ma_malloc(env, sizeof (Bittab*) * csa->number_of_sas);
  right = env_ma_malloc(env, sizeof (Bittab*) * csa->number_of_sas);
  L     = env_ma_malloc(env, sizeof (Bittab*) * csa->number_of_sas);
  R     = env_ma_malloc(env, sizeof (Bittab*) * csa->number_of_sas);

  for (i = 0; i < csa->number_of_sas; i++) {
    C[i]     = bittab_new(csa->number_of_sas, env);
    left[i]  = bittab_new(csa->number_of_sas, env);
    right[i] = bittab_new(csa->number_of_sas, env);
    L[i]     = bittab_new(csa->number_of_sas, env);
    R[i]     = bittab_new(csa->number_of_sas, env);
  }

  U_i      = bittab_new(csa->number_of_sas, env);
  SA_i     = bittab_new(csa->number_of_sas, env);
  SA_prime = bittab_new(csa->number_of_sas, env);

  splice_form = array_new(sizeof (unsigned long), env);

  /* compute sets */
  compute_C(C, csa, env);
  compute_left(left, csa, env);
  compute_right(right, csa, env);
  compute_L(L, C, left, csa->number_of_sas, env);
  compute_R(R, C, right, csa->number_of_sas, env);

  /* U_0 = SA */
  for (i = 0; i < csa->number_of_sas; i++)
    bittab_set_bit(U_i, i);

#ifndef NDEBUG
  /* preparation for assertion below */
  u_i_minus_1_size = bittab_count_set_bits(U_i);
#endif
  while (bittab_is_true(U_i)) {
    sa_i = UNDEF_ULONG;
    for (sa_prime  = bittab_get_first_bitnum(U_i);
         sa_prime != bittab_get_last_bitnum(U_i);
         sa_prime  = bittab_get_next_bitnum(U_i, sa_prime)) {
      if (sa_i == UNDEF_ULONG) {
        sa_i = sa_prime;
        bittab_or(SA_i, L[sa_i], R[sa_i]);
        sa_i_size = bittab_count_set_bits(SA_i);
      }
      else {
        bittab_or(SA_prime, L[sa_prime], R[sa_prime]);
        sa_prime_size = bittab_count_set_bits(SA_prime);
        if (sa_prime_size > sa_i_size) {
          sa_i = sa_prime;
          sa_i_size = sa_prime_size;
          bittab_equal(SA_i, SA_prime);
        }
      }
    }

    /* make sure the computed splice form is maximal w.r.t. to compatibility */
    assert(splice_form_is_valid(SA_i, csa, env));

    /* process splice form */
    if (csa->process_splice_form) {
      array_reset(splice_form);
      bittab_get_all_bitnums(SA_i, splice_form, env);
      csa->process_splice_form(splice_form, csa->set_of_sas, csa->number_of_sas,
                               csa->size_of_sa, csa->userdata, env);
    }

    /* U_i = U_i-1 \ SA_i */
    bittab_nand(U_i, U_i, SA_i);

#ifndef NDEBUG
    /* ensure that |U_i| < |U_i-1| */
    u_i_size = bittab_count_set_bits(U_i);
    assert(u_i_size < u_i_minus_1_size);
    u_i_minus_1_size = u_i_size;
#endif
  }

  /* free sets */
  for (i = 0; i < csa->number_of_sas; i++) {
    bittab_delete(C[i], env);
    bittab_delete(left[i], env);
    bittab_delete(right[i], env);
    bittab_delete(L[i], env);
    bittab_delete(R[i], env);
  }
  env_ma_free(C, env);
  env_ma_free(left, env);
  env_ma_free(right, env);
  env_ma_free(L, env);
  env_ma_free(R, env);
  bittab_delete(U_i, env);
  bittab_delete(SA_i, env);
  bittab_delete(SA_prime, env);
  array_delete(splice_form, env);
}

void consensus_sa(const void *set_of_sas, unsigned long number_of_sas,
                  size_t size_of_sa, GetGenomicRangeFunc get_genomic_range,
                  GetStrandFunc get_strand, GetExonsFunc get_exons,
                  ProcessSpliceFormFunc process_splice_form, void *userdata,
                  Env *env)
{
  ConsensusSA csa;
  assert(set_of_sas && number_of_sas && size_of_sa);
  assert(get_genomic_range && get_strand && get_exons);
  assert(set_of_sas_is_sorted(set_of_sas, number_of_sas, size_of_sa,
                              get_genomic_range, env));
  env_log_log(env, "csa number_of_sas=%lu", number_of_sas);

  /* init */
  csa.set_of_sas          = set_of_sas;
  csa.number_of_sas       = number_of_sas;
  csa.size_of_sa          = size_of_sa;
  csa.get_genomic_range   = get_genomic_range;
  csa.get_strand          = get_strand;
  csa.get_exons           = get_exons;
  csa.process_splice_form = process_splice_form;
  csa.userdata            = userdata;

  /* computation */
  compute_csas(&csa, env);

  env_log_log(env, "csa finished");
}
