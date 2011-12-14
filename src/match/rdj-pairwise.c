/*
  Copyright (c) 2009-2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2009-2011 Center for Bioinformatics, University of Hamburg

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

#include "core/ma.h"
#include "core/bioseq.h"
#include "core/encseq.h"
#include "core/progressbar.h"
#include "core/unused_api.h"
#include "match/rdj-ovlfind-bf.h"
#include "match/rdj-ovlfind-dp.h"
#include "match/rdj-ovlfind-kmp.h"
#include "match/rdj-ovlfind-kmp.h"
#include "match/rdj-spmlist.h"
#include "match/rdj-pairwise.h"
/* unit test: */
#include "match/rdj-ensure-output.h"
#include "core/fasta.h"

struct Read
{
  char* seq;
  unsigned long len;
  unsigned long seqnum;
  bool direct;
  gt_kmp_t *pi;
};

struct Data
{
  GtSpmproc proc;
  GtSpmprocA proc_a;
  void* procdata;
  struct Read *u, *v;
  GtOvlfindMode mode;
};

#define GT_RDJ_PAIRWISE_INIT_STRUCT_DATA(DATA, PROC, PROC_A, PROCDATA, U, V, \
    MODE) \
  {\
    (DATA).proc = (PROC);\
    (DATA).proc_a = (PROC_A);\
    (DATA).procdata = (PROCDATA);\
    (DATA).u = (U);\
    (DATA).v = (V);\
  }

static
void call_spmproc(unsigned long length, bool suffix_of_u, void *data)
{
  struct Data *d = data;
  if (suffix_of_u)
    d->proc(d->u->seqnum, d->v->seqnum, length,
        d->u->direct, d->v->direct, d->procdata);
  else
    d->proc(d->v->seqnum, d->u->seqnum, length, d->v->direct,
            d->u->direct, d->procdata);
}

static inline GtContfind find_exact_overlaps(struct Data *d, bool use_kmp,
    unsigned long min_length, bool find_nonmaximal)
{
  GtContfind retval;

  if (d->u->seqnum == d->v->seqnum && d->v->direct)
  {
    if (!(d->mode == GT_OVLFIND_SPM || d->mode == GT_OVLFIND_ALL))
      return GT_CONTFIND_EQ;
    if (use_kmp)
      retval = gt_ovlfind_kmp(d->u->seq, d->u->len, d->u->pi, NULL, 0, NULL,
          d->mode, min_length, find_nonmaximal, call_spmproc, d);
    else
      retval = gt_ovlfind_bf(d->u->seq, d->u->len, NULL, 0, d->mode,
          min_length, find_nonmaximal, call_spmproc, d);
  }
  else
  {
    if (use_kmp)
      retval = gt_ovlfind_kmp(d->u->seq, d->u->len, d->u->pi, d->v->seq,
          d->v->len, d->v->pi, d->mode, min_length, find_nonmaximal,
          call_spmproc, d);
    else
      retval = gt_ovlfind_bf(d->u->seq, d->u->len, d->v->seq, d->v->len,
          d->mode, min_length, find_nonmaximal, call_spmproc, d);
  }
  return retval;
}

static
void call_spmproc_a(unsigned long length_on_u, unsigned long length_on_v,
                    unsigned long unit_edist, bool suffix_of_u, void *data)
{
  struct Data *d = data;
  if (suffix_of_u)
    d->proc_a(d->u->seqnum, d->v->seqnum, length_on_u, length_on_v,
              unit_edist, d->u->direct, d->v->direct, d->procdata);
  else
    d->proc_a(d->v->seqnum, d->u->seqnum, length_on_v, length_on_u,
              unit_edist, d->v->direct, d->u->direct, d->procdata);
}

static inline
GtContfind find_approx_overlaps(struct Data *d, double max_error,
                                unsigned long min_length, bool find_nonmaximal)
{
  GtContfind retval;

  if (d->u->seqnum == d->v->seqnum && d->v->direct)
  {
    if (!(d->mode == GT_OVLFIND_SPM || d->mode == GT_OVLFIND_ALL))
      return GT_CONTFIND_EQ;
    retval = gt_ovlfind_dp(d->u->seq, d->u->len, NULL, 0, max_error, d->mode,
        min_length, find_nonmaximal, call_spmproc_a, d);
  }
  else
  {
    retval = gt_ovlfind_dp(d->u->seq, d->u->len, d->v->seq, d->v->len,
        max_error, d->mode, min_length, find_nonmaximal, call_spmproc_a, d);
  }
  return retval;
}

static
gt_kmp_t** prepare_kmp_values(const GtEncseq *encseq, const GtBioseq *bs,
    unsigned long n)
{
  unsigned long i, len, startpos;
  char *seq;
  gt_kmp_t **kmp_values;

  kmp_values = gt_malloc(sizeof (gt_kmp_t*) * n);
  gt_assert(encseq != NULL || bs != NULL);
  if (encseq == NULL)
  {
    gt_assert(bs != NULL);
    for (i = 0; i < n; i++)
      kmp_values[i] = gt_kmp_preproc(
                        gt_seq_get_orig(gt_bioseq_get_seq((GtBioseq*)bs, i)),
                        gt_bioseq_get_sequence_length(bs, i));
  }
  else
  {
    for (i = 0; i < n; i++)
    {
      len = gt_encseq_seqlength(encseq, i);
      seq = gt_malloc(sizeof (char) * (len + 1));
      startpos = gt_encseq_seqstartpos(encseq, i);
      gt_encseq_extract_decoded(encseq, seq, startpos, startpos + len - 1);
      seq[len] = '\0';

      kmp_values[i] = gt_kmp_preproc(seq, len);
      gt_free(seq);
    }
  }
  return kmp_values;
}

static
void free_kmp_values(gt_kmp_t** kmp_values, unsigned long n)
{
  unsigned long i;
  gt_assert(kmp_values != NULL);
  for (i = 0; i < n; i++) gt_free(kmp_values[i]);
  gt_free(kmp_values);
}

static inline
void mark_contained(GtContfind c, unsigned long u_seqnum,
                    unsigned long v_seqnum, GtBitsequence *cntreads)
{
  gt_assert(cntreads != NULL);
  switch (c)
  {
    case GT_CONTFIND_U:
      GT_SETIBIT(cntreads, u_seqnum);
      break;
    case GT_CONTFIND_EQ:
      /* convention: skip one with higher seqnum (which is v) */
      if (v_seqnum != u_seqnum)
        GT_SETIBIT(cntreads, v_seqnum);
      break;
    case GT_CONTFIND_V:
      GT_SETIBIT(cntreads, v_seqnum);
      break;
    default:
      /* nothing to do */
      break;
  }
}

static inline
int rdj_pairwise_generic(bool use_dp, GtOvlfindMode m, GtEncseq *encseq,
                           const char *filename, bool revcompl,
                           bool show_progressbar, bool use_kmp,
                           double max_error, unsigned long min_length,
                           bool find_nonmaximal, GtSpmproc proc,
                           GtSpmprocA proc_a, void* procdata,
                           bool cntfilter, GtBitsequence *cntreads_in,
                           GtBitsequence **cntreads_out,
                           unsigned long *nofreads, GtError *err)
{
  GtContfind containment_status;
  GtBitsequence *cntreads = NULL;
  GtBioseq *bs = NULL;
  unsigned long long progress = 0;
  unsigned long i, j, startpos, v_seqnum, nofsequences, n;
  struct Read u, v;
  struct Data d;
  gt_kmp_t** kmp_values = NULL;

  GT_RDJ_PAIRWISE_INIT_STRUCT_DATA(d, proc, proc_a, procdata, &u, &v, 0);

  /* encseq XOR filename required: */
  gt_assert((encseq != NULL) || (filename != NULL));
  gt_assert(!((encseq != NULL) && (filename != NULL)));

  d.mode = m;
  if ((m == GT_OVLFIND_ALL) && cntfilter)
    d.mode = GT_OVLFIND_PROPER_SPM;

  if (encseq == NULL)
  {
    bs = gt_bioseq_new(filename, err);
    if (bs == NULL) return 1;
    n = gt_bioseq_number_of_sequences(bs);
    if (use_kmp)
      kmp_values = prepare_kmp_values(NULL, bs, n);
  }
  else
  {
    n = gt_encseq_num_of_sequences(encseq);
    if (use_kmp)
      kmp_values = prepare_kmp_values(encseq, NULL, n);
  }
  nofsequences = n;
  if (revcompl)
    n = n >> 1;
  if (cntreads_in != NULL)
    cntreads = cntreads_in;
  else if (m != GT_OVLFIND_SPM)
    GT_INITBITTAB(cntreads, n);
  if (show_progressbar) gt_progressbar_start(&progress, (unsigned long long)n *
      ((unsigned long long)n - 1ULL) / 2ULL);

  for (i = 0; i < n; i++)
  {
    u.seqnum = i;
    u.direct = true;
    if (encseq == NULL)
    {
      u.len = gt_bioseq_get_sequence_length(bs, i);
      u.seq = (char*) gt_seq_get_orig(gt_bioseq_get_seq(bs, i));
    }
    else
    {
      u.len = gt_encseq_seqlength(encseq, i);
      u.seq = gt_malloc(sizeof (char) * (u.len + 1));
      startpos = gt_encseq_seqstartpos(encseq, i);
      gt_encseq_extract_decoded(encseq, u.seq, startpos, startpos + u.len - 1);
      u.seq[u.len] = '\0';
    }
    if (use_kmp)
    {
      gt_assert(kmp_values != NULL);
      u.pi = kmp_values[i];
    }

    for (j = i; j < n; j++)
    {
      if (cntfilter)
      {
        gt_assert(cntreads != NULL);
        if ((bool)GT_ISIBITSET(cntreads, i)) break;
        if ((bool)GT_ISIBITSET(cntreads, j)) continue;
      }

      v.seqnum = j;

      /* find overlaps using direct v */
      v.direct = true;
      if (encseq == NULL)
      {
        v.len = gt_bioseq_get_sequence_length(bs, j);
        v.seq = (char*) gt_seq_get_orig(gt_bioseq_get_seq(bs, j));
      }
      else
      {
        v.len = gt_encseq_seqlength(encseq, j);
        v.seq = gt_malloc(sizeof (char) * (v.len + 1));
        startpos = gt_encseq_seqstartpos(encseq, j);
        gt_encseq_extract_decoded(encseq, v.seq, startpos,
          startpos + v.len - 1);
        v.seq[v.len] = '\0';
      }
      if (use_kmp)
      {
        gt_assert(kmp_values != NULL);
        v.pi = kmp_values[j];
      }
      containment_status = use_dp
          ? find_approx_overlaps(&d, max_error, min_length, find_nonmaximal)
          : find_exact_overlaps(&d, use_kmp, min_length, find_nonmaximal);
      if (m != GT_OVLFIND_SPM)
        mark_contained(containment_status, u.seqnum, v.seqnum, cntreads);

      /* find overlaps using reverse complement of v */
      if (revcompl)
      {
        v_seqnum =  nofsequences - j - 1;
        v.direct = false;
        if (encseq == NULL)
        {
          v.seq = (char*) gt_seq_get_orig(gt_bioseq_get_seq(bs, v_seqnum));
        }
        else
        {
          gt_assert(gt_encseq_seqlength(encseq, j) ==
            gt_encseq_seqlength(encseq, v_seqnum));
          startpos = gt_encseq_seqstartpos(encseq, v_seqnum);
          gt_encseq_extract_decoded(encseq, v.seq, startpos,
            startpos + v.len - 1);
        }
        if (use_kmp)
        {
          gt_assert(kmp_values != NULL);
          v.pi = kmp_values[v_seqnum];
        }
        containment_status = use_dp
          ? find_approx_overlaps(&d, max_error, min_length, find_nonmaximal)
          : find_exact_overlaps(&d, use_kmp, min_length, find_nonmaximal);
        if (m != GT_OVLFIND_SPM)
          mark_contained(containment_status, u.seqnum, v.seqnum, cntreads);
      }
      if (encseq != NULL)
        gt_free(v.seq);
      progress++;
    }
    if (encseq != NULL)
      gt_free(u.seq);
  }

  if (cntreads_out != NULL)
    *cntreads_out = cntreads;
  else if (cntreads_in == NULL)
    gt_free(cntreads);
  if (nofreads != NULL)
    *nofreads = n;
  gt_bioseq_delete(bs);
  if (use_kmp)
    free_kmp_values(kmp_values, revcompl ? n << 1 : n);
  if (show_progressbar)
    gt_progressbar_stop();
  return 0;
}

#ifndef NDEBUG
static inline
bool rdj_pairwise_check_arguments(GtOvlfindMode m, void* proc, void* procdata,
                     GtBitsequence *cntreads_in, GtBitsequence **cntreads_out,
                     bool cntfilter, GtError *err)
{
  /*@i1@*/ gt_error_check(err);
  gt_assert(m == GT_OVLFIND_CNT || m == GT_OVLFIND_ALL || m == GT_OVLFIND_SPM);
  if (proc == NULL)
  {
    gt_assert(m == GT_OVLFIND_CNT);
    gt_assert(procdata == NULL);
  }
  else
    gt_assert(m == GT_OVLFIND_SPM || m == GT_OVLFIND_ALL);
  if (cntfilter || (cntreads_out != NULL))
    gt_assert((m == GT_OVLFIND_CNT) ||
              (m == GT_OVLFIND_ALL) ||
              (cntreads_in != NULL));
  return true; /* allow to wrap in gt_assert() */
}
#endif

int gt_rdj_pairwise_exact(GtOvlfindMode m, GtEncseq *encseq,
                            const char *filename, bool revcompl,
                            bool show_progressbar, bool use_kmp,
                            unsigned long min_length, bool find_nonmaximal,
                            GtSpmproc proc, void* procdata,
                            bool cntfilter, GtBitsequence *cntreads_in,
                            GtBitsequence **cntreads_out,
                            unsigned long *nofreads, GtError *err)
{
  gt_assert(rdj_pairwise_check_arguments(m, proc, procdata, cntreads_in,
                            cntreads_out, cntfilter, err));
  return rdj_pairwise_generic(false, m, encseq, filename, revcompl,
    show_progressbar, use_kmp, 0.0, min_length, find_nonmaximal, proc, NULL,
    procdata, cntfilter, cntreads_in, cntreads_out, nofreads, err);
}

int gt_rdj_pairwise_approx(GtOvlfindMode m,  GtEncseq *encseq,
                             const char *filename, bool revcompl,
                             bool show_progressbar, double max_error,
                             unsigned long min_length, bool find_nonmaximal,
                             GtSpmprocA proc, void* procdata,
                             bool cntfilter, GtBitsequence *cntreads_in,
                             GtBitsequence **cntreads_out,
                             unsigned long *nofreads, GtError *err)
{
  gt_assert(rdj_pairwise_check_arguments(m, proc, procdata, cntreads_in,
                            cntreads_out, cntfilter, err));
  return rdj_pairwise_generic(true, m, encseq, filename, revcompl,
    show_progressbar, false, max_error, min_length, find_nonmaximal, NULL,
    proc, procdata, cntfilter, cntreads_in, cntreads_out, nofreads, err);
}

int gt_rdj_pairwise_unit_test(GtError *err)
{
  int had_err = 0;
  int retval;
  GT_ENSURE_OUTPUT_DECLARE(200);
  GtStr *fas_str, *bs_str;
  GtFile *fas_file;

  /*@i1@*/ gt_error_check(err);

  /* prepare fasta file with some reads */
  fas_str = gt_str_new();
  fas_file = gt_file_new_from_fileptr(gt_xtmpfp(fas_str));
  gt_fasta_show_entry("0", "caacgaatctgc", 12UL, 60UL, fas_file);
  gt_fasta_show_entry("1", "cacacccacaac", 12UL, 60UL, fas_file);
  gt_fasta_show_entry("2", "tgcacgtggcca", 12UL, 60UL, fas_file);
  gt_fasta_show_entry("3", "caacgcagctgg", 12UL, 60UL, fas_file);
  gt_fasta_show_entry("4", "gtgcactagtga", 12UL, 60UL, fas_file);
  gt_fasta_show_entry("4rc", "tcactagtgcac", 12UL, 60UL, fas_file);
  gt_fasta_show_entry("3rc", "ccagctgcgttg", 12UL, 60UL, fas_file);
  gt_fasta_show_entry("2rc", "tggccacgtgca", 12UL, 60UL, fas_file);
  gt_fasta_show_entry("1rc", "gttgtgggtgtg", 12UL, 60UL, fas_file);
  gt_fasta_show_entry("0rc", "gcagattcgttg", 12UL, 60UL, fas_file);
  gt_file_delete(fas_file);

  /* run rdj_pairwise */
  GT_ENSURE_OUTPUT(
      retval = gt_rdj_pairwise_exact(GT_OVLFIND_SPM, NULL, gt_str_get(fas_str),
        true, false, false, 3UL, false, gt_spmproc_show_ascii, outfp, false,
        NULL, NULL, NULL, err),
      "1 + 0 + 4\n0 + 2 + 3\n1 + 3 + 4\n4 - 1 + 3\n"
      "2 + 2 - 6\n2 - 2 + 4\n2 + 3 - 3\n4 - 2 + 5\n");
  gt_ensure(had_err, retval == 0);

  GT_ENSURE_OUTPUT(
      retval = gt_rdj_pairwise_approx(GT_OVLFIND_SPM, NULL, gt_str_get(fas_str),
        true, false, 0.0, 3UL, false, gt_spmproc_a_show_ascii, outfp,
        false, NULL, NULL, NULL, err),
      "1 + 0 + 4 4 0\n0 + 2 + 3 3 0\n1 + 3 + 4 4 0\n4 - 1 + 3 3 0\n"
      "2 + 2 - 6 6 0\n2 - 2 + 4 4 0\n2 + 3 - 3 3 0\n4 - 2 + 5 5 0\n");
  gt_ensure(had_err, retval == 0);

  /* clean fasta and bioseq files */
  bs_str = gt_str_clone(fas_str);
  gt_str_append_cstr(bs_str, GT_BIOSEQ_INDEX);
  gt_xremove(gt_str_get(bs_str));
  gt_str_delete(bs_str);
  bs_str = gt_str_clone(fas_str);
  gt_str_append_cstr(bs_str, GT_BIOSEQ_RAW);
  gt_xremove(gt_str_get(bs_str));
  gt_str_delete(bs_str);
  gt_xremove(gt_str_get(fas_str));
  gt_str_delete(fas_str);

  if (had_err == 0)
    had_err = gt_ovlfind_bf_unit_test(err);
  if (had_err == 0)
    had_err = gt_ovlfind_kmp_unit_test(err);
  if (had_err == 0)
    had_err = gt_ovlfind_dp_unit_test(err);

  return had_err;
}
