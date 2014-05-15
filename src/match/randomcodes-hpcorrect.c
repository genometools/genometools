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

#include "core/fa.h"
#include "core/log_api.h"
#include "core/mathsupport.h"
#include "core/xansi_api.h"
#include "core/undef_api.h"
#include "extended/union_find.h"
#include "match/randomcodes-hpcorrect.h"
#include "match/triangular_def.h"

//#define GT_RANDOMCODES_HPCORRECT_DEBUG
//#define GT_RANDOMCODES_HPCORRECT_VERBOSE

#define GT_RANDOMCODES_HPCORRECT_HMERS_INIT  128UL
#define GT_RANDOMCODES_HPCORRECT_HMERS_EXTRA 128UL

#define GT_RANDOMCODES_HPCORRECT_PW_SCORE(SDATA,I,J) \
  ((SDATA)->pw_scores[\
     GT_STRICT_TRIANGULAR_SADDR((SDATA)->hmers_alloc,(I),(J))])

/* ptr to i-th hmer */
#define GT_RANDOMCODES_HPCORRECT_HMER(SDATA, I) \
  ((SDATA)->hmers + ((I) * (SDATA)->hmers_width))

/* define a type for hmer pw_scores */
typedef GtWord gt_rchc_score_t;
#ifndef _WIN64
#define GT_RCHC_PRIscore "ld"
#else
#define GT_RCHC_PRIscore "I64d"
#endif

static int gt_rchc_score_t_cmp(const void *a, const void *b)
{
  const gt_rchc_score_t *sa = a;
  const gt_rchc_score_t *sb = b;
  return (int)((*sa > *sb) - (*sa < *sb));
}

#define GT_RANDOMCODES_HPCORRECT_CLUSTERNUM(SDATA,N) \
  ((SDATA)->union_find != NULL \
    ? gt_union_find_find((SDATA)->union_find, (N)) \
    : (SDATA)->hmer_cluster[N])

struct GtRandomcodesHpcorrectData {
  const GtEncseq *encseq;
  GtHplstore *hplstore;
  GtUword nofreads;
  GtUword mirror_nofseqs;
  GtUword k;
  unsigned int maxuntrustp;
  unsigned int mintrustp;
  uint8_t *hmers;
  GtUword *cluster_size;
  GtUword hmers_alloc;
  GtUword hmers_width;
  uint8_t *consensus_hmer;
  unsigned int *hdistri;
  bool skip_read_ends;
  GtUword frompos;
  GtUword topos;
  bool skip_rc;
  GtWord clustering_minscore;
  GtUword *hmer_cluster;
  bool *skip;
  GtUnionFind *union_find;
  GtFile *outfile;
  gt_rchc_score_t *pw_scores;
  gt_rchc_score_t *pw_scores_copy;
  unsigned int quantile;
  GtUword maxwidth;
  bool manhattan;
  /* firstpass data: */
  bool firstpass;
  unsigned int r;
  GtUword *rextset;
  unsigned int *rextset_size;
  gt_rchc_score_t rext_cl_minscore, rext_I_minscore, rext_J_minscore,
                  rext_R_minscore, rext_D_minscore, rext_J_lminscore;
  GtUword rext_J_lwidth;
  GtEncseqReader *esr1, *esr2;
};

#define GT_RCHC_HSCORE_MAXH (uint8_t)9
#define GT_RCHC_HSCORE_BEST (gt_rchc_score_t)4L
static gt_rchc_score_t gt_rchc_hscore[100] =
   { 0L, -2L, -4L, -4L, -5L, -6L, -6L, -6L, -6L, -6L,
    -2L,  1L, -1L, -3L, -4L, -4L, -4L, -6L, -6L, -6L,
    -4L, -1L,  2L,  0L, -2L, -3L, -3L, -6L, -6L, -6L,
    -4L, -3L,  0L,  2L,  1L, -1L, -1L, -1L, -2L, -2L,
    -5L, -4L, -2L,  1L,  2L,  1L,  0L,  0L,  0L,  0L,
    -6L, -4L, -3L, -1L,  1L,  3L,  2L,  1L,  1L,  1L,
    -6L, -4L, -3L, -1L,  0L,  2L,  3L,  3L,  2L,  2L,
    -6L, -6L, -6L, -1L,  0L,  1L,  3L,  4L,  3L,  3L,
    -6L, -6L, -6L, -2L,  0L,  1L,  2L,  3L,  4L,  3L,
    -6L, -6L, -6L, -2L,  0L,  1L,  2L,  3L,  3L,  4L};

static inline gt_rchc_score_t gt_randomcodes_hpcorrect_compute_hscore_manhattan(
    bool *identical, const uint8_t *hmer1, const uint8_t *hmer2,
    GtUword k)
{
  GtUword i;
  gt_rchc_score_t total_score = 0;
  if (identical != NULL)
    *identical = true;
  for (i = 0; i < k; i++)
  {
    uint8_t h1 = hmer1[i], h2 = hmer2[i];
    if (h1!=h2)
    {
      if (identical != NULL)
        *identical = false;
      total_score -= (h1 > h2 ? h1 - h2 : h2 - h1);
    }
  }
  return total_score;
}

#if 0
GT_UNUSED
static inline gt_rchc_score_t
    gt_randomcodes_hpcorrect_compute_hscore_manhattan_threshold(
    bool *identical, const uint8_t *hmer1, const uint8_t *hmer2,
    GtUword k, gt_rchc_score_t minscore /* use LONG_MAX to disable */)
{
  GtUword i;
  gt_rchc_score_t total_score = 0;
  gt_assert(identical != NULL);
  *identical = true;
  for (i = 0; i < k && (*identical || total_score >= minscore); i++)
  {
    uint8_t h1 = hmer1[i], h2 = hmer2[i];
    if (h1!=h2)
    {
      *identical = false;
      total_score -= (h1 > h2 ? h1 - h2 : h2 - h1);
    }
  }
  return total_score;
}

GT_UNUSED
static inline gt_rchc_score_t gt_randomcodes_hpcorrect_compute_hscore_flexible(
    bool *identical, const uint8_t *hmer1, const uint8_t *hmer2,
    GtUword k, bool manhattan,
    gt_rchc_score_t minscore /* use LONG_MAX to disable */)
{
  GtUword i;
  gt_rchc_score_t total_score = 0;
  gt_rchc_score_t max_reachable = manhattan ? 0 : (GT_RCHC_HSCORE_BEST * k);
  gt_assert(identical != NULL);
  *identical = true;
  for (i = 0; i < k && (*identical || max_reachable >= minscore); i++)
  {
    uint8_t h1 = hmer1[i], h2 = hmer2[i];
    if (h1!=h2)
      *identical = false;
    if (manhattan)
    {
      gt_rchc_score_t distance;
      distance = h1 > h2 ? h1 - h2 : h2 - h1;
      total_score -= distance;
      max_reachable -= distance;
    }
    else
    {
      gt_rchc_score_t score;
      if (h1 > GT_RCHC_HSCORE_MAXH)
        h1 = GT_RCHC_HSCORE_MAXH;
      if (h2 > GT_RCHC_HSCORE_MAXH)
        h2 = GT_RCHC_HSCORE_MAXH;
      score = gt_rchc_hscore[(GT_RCHC_HSCORE_MAXH+1)*h1+h2];
      total_score += score;
      max_reachable -= (GT_RCHC_HSCORE_BEST - score);
    }
  }
  return total_score;
}

GT_UNUSED
static inline gt_rchc_score_t
    gt_randomcodes_hpcorrect_compute_hscore_matrix_threshold(
    bool *identical, const uint8_t *hmer1, const uint8_t *hmer2,
    GtUword k, gt_rchc_score_t minscore /* use LONG_MAX to disable */)
{
  GtUword i;
  gt_rchc_score_t score, total_score = 0;
  gt_rchc_score_t max_reachable = GT_RCHC_HSCORE_BEST * k;
  gt_assert(identical != NULL);
  *identical = true;
  for (i = 0; i < k && (*identical || max_reachable >= minscore); i++)
  {
    uint8_t h1 = hmer1[i], h2 = hmer2[i];
    if (h1!=h2)
      *identical = false;
    if (h1 > GT_RCHC_HSCORE_MAXH)
      h1 = GT_RCHC_HSCORE_MAXH;
    if (h2 > GT_RCHC_HSCORE_MAXH)
      h2 = GT_RCHC_HSCORE_MAXH;
    score = gt_rchc_hscore[(GT_RCHC_HSCORE_MAXH+1)*h1+h2];
    total_score += score;
    max_reachable -= (GT_RCHC_HSCORE_BEST - score);
  }
  return total_score;
}

#endif
static inline gt_rchc_score_t
    gt_randomcodes_hpcorrect_compute_hscore(
    bool *identical, const uint8_t *hmer1, const uint8_t *hmer2,
    GtUword k)
{
  GtUword i;
  gt_rchc_score_t total_score = 0;
  gt_assert(identical != NULL);
  *identical = true;
  for (i = 0; i < k; i++)
  {
    uint8_t h1 = hmer1[i], h2 = hmer2[i];
    if (h1!=h2)
      *identical = false;
    if (h1 > GT_RCHC_HSCORE_MAXH)
      h1 = GT_RCHC_HSCORE_MAXH;
    if (h2 > GT_RCHC_HSCORE_MAXH)
      h2 = GT_RCHC_HSCORE_MAXH;
    total_score += gt_rchc_hscore[(GT_RCHC_HSCORE_MAXH+1)*h1+h2];
  }
  return total_score;
}

GT_UNUSED
static void gt_randomcodes_hpcorrect_show_expanded_kmer(
    GtRandomcodesHpcorrectData *sdata, const GtUword *bucketofsuffixes,
    const GtSeqnumrelpos *snrp, GtUword sfx, GtUword info)
{
  GtUword relpos, seqnum, startpos, kmerpos;
  relpos = gt_seqnumrelpos_decode_relpos(snrp, bucketofsuffixes[sfx]);
  seqnum = gt_seqnumrelpos_decode_seqnum(snrp, bucketofsuffixes[sfx]);
  startpos = gt_encseq_seqstartpos(sdata->encseq, seqnum);
  kmerpos = startpos + relpos;
  gt_file_xprintf(sdata->outfile, "# ["GT_WU"", sfx);
  if (info != GT_UNDEF_UWORD)
    gt_file_xprintf(sdata->outfile, "-"GT_WU"", info);
  gt_file_xprintf(sdata->outfile, "] ");
  gt_hplstore_show_decoded_sequence(sdata->outfile, sdata->hplstore,
      sdata->encseq, kmerpos, sdata->k);
  if (sdata->firstpass)
  {
    gt_file_xprintf(sdata->outfile, "|");
    gt_hplstore_show_decoded_sequence(sdata->outfile, sdata->hplstore,
        sdata->encseq, kmerpos + sdata->k, (GtUword)sdata->r + 1UL);
  }
  gt_file_xprintf(sdata->outfile, " ("GT_WU":"GT_WU")\n", seqnum, relpos);
}

GT_UNUSED
static inline void gt_randomcodes_hpcorrect_show_all_expanded_kmers(
    const GtSeqnumrelpos *snrp, const GtUword *suffixes,
    GtUword nofsuffixes, GtRandomcodesHpcorrectData *sdata,
    GtUword clusternum /* GT_UNDEF_UWORD for all clusters */)
{
  GtUword i;
  for (i = 0; i < nofsuffixes; i++)
  {
    if (clusternum == GT_UNDEF_UWORD || \
        GT_RANDOMCODES_HPCORRECT_CLUSTERNUM(sdata, i) == clusternum)
      gt_randomcodes_hpcorrect_show_expanded_kmer(sdata, suffixes, snrp, i,
          clusternum);
  }
}

static inline uint8_t gt_randomcodes_hpcorrect_compute_consensus_for_pos(
    GtUword pos, GtUword nofsuffixes,
    GtRandomcodesHpcorrectData *sdata, GtUword clusternum)
{
  GtUword i;
  unsigned int maxfreq = 0;
  uint8_t value, consensus = GT_HPLSTORE_UNDEF;
  unsigned int mintrust = (unsigned int)sdata->cluster_size[clusternum] *
    sdata->mintrustp / 100U;
  unsigned int *hdistri = sdata->hdistri + pos * GT_HPLSTORE_RANGE;
  for (i = 0; (uint8_t)i < GT_HPLSTORE_RANGE; i++)
    hdistri[i] = 0;
  for (i = 0; i < nofsuffixes; i++)
  {
    if (GT_RANDOMCODES_HPCORRECT_CLUSTERNUM(sdata, i) == clusternum)
    {
      value = GT_RANDOMCODES_HPCORRECT_HMER(sdata, i)[pos];
      hdistri[value]++;
      if (hdistri[value] > maxfreq)
      {
        consensus = value;
        maxfreq = hdistri[value];
      }
      else if (hdistri[value] == maxfreq)
      {
        consensus = GT_HPLSTORE_UNDEF;
      }
    }
  }
  if (maxfreq < mintrust)
    consensus = GT_HPLSTORE_UNDEF;
  return consensus;
}

static inline void gt_randomcodes_hpcorrect_compute_consensus_hmer(
    GtUword nofsuffixes, GtRandomcodesHpcorrectData *sdata,
    GtUword clusternum)
{
  GtUword pos;
  for (pos = 0; pos < sdata->k; pos++)
  {
    sdata->consensus_hmer[pos] =
      gt_randomcodes_hpcorrect_compute_consensus_for_pos(pos, nofsuffixes,
          sdata, clusternum);
  }
}

static inline void gt_randomcodes_hpcorrect_fill_pw_scores_matrix(
    GtUword nofsuffixes, GtRandomcodesHpcorrectData *sdata)
{
  GtUword i, j;
  bool identical;
  for (i = 0; i < nofsuffixes; i++)
  {
    for (j = i + 1; j < nofsuffixes; j++)
    {
      GT_RANDOMCODES_HPCORRECT_PW_SCORE(sdata,i,j) =
        gt_randomcodes_hpcorrect_compute_hscore(
            &identical, GT_RANDOMCODES_HPCORRECT_HMER(sdata,i),
            GT_RANDOMCODES_HPCORRECT_HMER(sdata,j), sdata->k
            /*, sdata->manhattan, LONG_MAX*/);
    }
  }
}

GT_UNUSED
static inline void gt_randomcodes_hpcorrect_show_pw_scores(
    GtUword nofsuffixes, GtRandomcodesHpcorrectData *sdata)
{
  GtUword i, j;
  gt_rchc_score_t score;
  gt_file_xprintf(sdata->outfile, "# ");
  for (i = 0; i < nofsuffixes; i++)
  {
    gt_file_xprintf(sdata->outfile, "\t("GT_WU")", i);
  }
  gt_file_xprintf(sdata->outfile, "\n");
  for (i = 0; i < nofsuffixes; i++)
  {
    gt_file_xprintf(sdata->outfile, "# ("GT_WU")\t", i);
    for (j = 0; j <= i; j++)
    {
      gt_file_xprintf(sdata->outfile, "\t");
    }
    for (j = i + 1; j < nofsuffixes; j++)
    {
      score = GT_RANDOMCODES_HPCORRECT_PW_SCORE(sdata,i,j);
      gt_file_xprintf(sdata->outfile, "%3"GT_RCHC_PRIscore"\t", score);
    }
    gt_file_xprintf(sdata->outfile, "\n");
  }
}

static inline GtUword gt_randomcodes_hpcorrect_cluster_greedy(
    bool *allidentical, GtUword nofsuffixes, GtWord clustering_minscore,
    GtRandomcodesHpcorrectData *sdata)
{
  gt_rchc_score_t hscore;
  GtUword i, j, nofclusters = 0;
  bool identical = false;
  gt_assert(allidentical != NULL);
  *allidentical = true;
  gt_assert(sdata->hmer_cluster != NULL);
  for (i = 0; i < nofsuffixes; i++)
    sdata->hmer_cluster[i] = GT_UNDEF_UWORD;
  (void)memset(sdata->cluster_size, 0, sizeof (sdata->cluster_size) *
      nofsuffixes);
  for (i = 0; i < nofsuffixes; i++)
  {
    if (sdata->hmer_cluster[i] == GT_UNDEF_UWORD)
      sdata->hmer_cluster[i] = (++nofclusters);
    for (j = i+1; j < nofsuffixes; j++)
    {
      hscore = (sdata->pw_scores == NULL) ?
        gt_randomcodes_hpcorrect_compute_hscore(
            &identical, GT_RANDOMCODES_HPCORRECT_HMER(sdata,i),
            GT_RANDOMCODES_HPCORRECT_HMER(sdata,j), sdata->k
            /*, sdata->manhattan, clustering_minscore*/)
        : GT_RANDOMCODES_HPCORRECT_PW_SCORE(sdata,i,j);
      if (!identical)
        *allidentical = false;
      if (hscore >= clustering_minscore)
      {
        gt_assert(sdata->hmer_cluster[i] != GT_UNDEF_UWORD);
        sdata->hmer_cluster[j] = sdata->hmer_cluster[i];
      }
    }
    sdata->cluster_size[sdata->hmer_cluster[i]]++;
  }
  return nofclusters;
}

static inline GtUword gt_randomcodes_hpcorrect_cluster_union_find(
    bool *allidentical, GtUword nofsuffixes, GtWord clustering_minscore,
    GtRandomcodesHpcorrectData *sdata)
{
  gt_rchc_score_t hscore;
  GtUword i, j, nofclusters;
  bool identical = false;
  gt_assert(allidentical != NULL);
  *allidentical = true;
  gt_assert(sdata->union_find != NULL);
  gt_union_find_reset(sdata->union_find, nofsuffixes);
  gt_assert(sdata->skip != NULL);
  (void)memset(sdata->skip, 0, sizeof (sdata->skip) * nofsuffixes);
  for (i = 0; i < nofsuffixes; i++)
  {
    if (sdata->skip[i])
      continue;
    for (j = i+1; j < nofsuffixes; j++)
    {
      if (sdata->skip[j])
        continue;
      hscore = (sdata->pw_scores == NULL) ?
        gt_randomcodes_hpcorrect_compute_hscore(
            &identical, GT_RANDOMCODES_HPCORRECT_HMER(sdata,i),
            GT_RANDOMCODES_HPCORRECT_HMER(sdata,j), sdata->k
            /*,sdata->manhattan, clustering_minscore*/)
        : GT_RANDOMCODES_HPCORRECT_PW_SCORE(sdata,i,j);
      if (!identical)
        *allidentical = false;
      if (hscore >= clustering_minscore)
      {
        gt_union_find_union(sdata->union_find, i, j);
        if (identical)
          sdata->skip[j] = true;
      }
    }
  }
  (void)memset(sdata->cluster_size, 0, sizeof (sdata->cluster_size) *
      nofsuffixes);
  for (i = 0; i < nofsuffixes; i++)
    sdata->cluster_size[gt_union_find_find(sdata->union_find, i)]++;
  nofclusters = 0;
  for (i = 0; i < nofsuffixes; i++)
    if (sdata->cluster_size[i] > 0)
      nofclusters++;
  return nofclusters;
}

GT_UNUSED
static inline void gt_randomcodes_hpcorrect_show_kplus1(
    const GtSeqnumrelpos *snrp, const GtUword *suffixes,
    GtUword nofsuffixes, GtRandomcodesHpcorrectData *sdata,
    GtUword clusternum)
{
  GtUword i, relpos, seqnum, startpos, kmerpos, kplus1pos;
  char kplus1char;
  bool rc;
  gt_file_xprintf(sdata->outfile, "# clusternum: "GT_WU"\n", clusternum);
  for (i = 0; i < nofsuffixes; i++)
  {
    if (GT_RANDOMCODES_HPCORRECT_CLUSTERNUM(sdata, i) != clusternum)
      continue;
    relpos = gt_seqnumrelpos_decode_relpos(snrp, suffixes[i]);
    seqnum = gt_seqnumrelpos_decode_seqnum(snrp, suffixes[i]);
    startpos = gt_encseq_seqstartpos(sdata->encseq, seqnum);
    kmerpos = startpos + relpos;
    kplus1char = gt_encseq_get_decoded_char(sdata->encseq, kmerpos + sdata->k,
        GT_READMODE_FORWARD);
    kplus1pos = relpos + sdata->k;
    rc = (seqnum >= sdata->nofreads);
    if (rc)
    {
      GtUword seqlen;
      if (sdata->skip_rc)
        continue;
      seqlen = gt_encseq_seqlength(sdata->encseq, seqnum);
      gt_file_xprintf(sdata->outfile,"# rc: "GT_WU" "GT_WU", seqlen="GT_WU"\n",
          seqnum, kplus1pos, seqlen);
      kplus1pos = seqlen - 1UL - kplus1pos;
      seqnum = sdata->mirror_nofseqs - 1UL - seqnum;
    }
    gt_file_xprintf(sdata->outfile,"# kplus1\t"GT_WU"\t"GT_WU"\t%c\trextset="
                    GT_WU"\n", seqnum, kplus1pos, kplus1char,
                    sdata->rextset[i]);
  }
}

static char gt_randomcodes_complement(char dna_char)
{
  switch (dna_char) {
    case 'a': case 'A': return 'T';
    case 'c': case 'C': return 'G';
    case 'g': case 'G': return 'C';
    case 't': case 'T': return 'A';
  }
  return 'n';
}

/* compare sdata->r letters from position kmerpos + sdata->k + offset */
static inline bool gt_randomcodes_hpcorrect_same_rext_letters(
    GtRandomcodesHpcorrectData *sdata, GtUword u_kmerpos, int u_offset,
    GtUword t_kmerpos, int t_offset)
{
  GtCommonunits commonunits;
  return (u_kmerpos + u_offset == t_kmerpos + t_offset) ||
      (gt_encseq_compare_viatwobitencoding(&commonunits,
        sdata->encseq, sdata->encseq, GT_READMODE_FORWARD, sdata->esr1,
        sdata->esr2, u_kmerpos + u_offset + sdata->k,
        t_kmerpos + t_offset + sdata->k, 0, (GtUword)sdata->r) == 0);
}

/* compare right extensions letters and homopolymer lengths */
static inline bool gt_randomcodes_hpcorrect_similar_rext(
    GtRandomcodesHpcorrectData *sdata,
    GtUword u_sfx, GtUword u_kmerpos, int u_offset,
    GtUword t_sfx, GtUword t_kmerpos, int t_offset,
    gt_rchc_score_t minscore)
{
  return (gt_randomcodes_hpcorrect_same_rext_letters(sdata, u_kmerpos, u_offset,
        t_kmerpos, t_offset)) &&
    (gt_randomcodes_hpcorrect_compute_hscore_manhattan(NULL,
          GT_RANDOMCODES_HPCORRECT_HMER(sdata, u_sfx) + sdata->k + u_offset,
          GT_RANDOMCODES_HPCORRECT_HMER(sdata, t_sfx) + sdata->k + t_offset,
          (GtUword)sdata->r) >= minscore);
}

static inline bool gt_randomcodes_hpcorrect_detect_double_insertion(
    GtRandomcodesHpcorrectData *sdata, GtUword u_sfx,
    GtUword u_kmerpos, GtUword t_sfx, GtUword t_kmerpos,
    gt_rchc_score_t r_minscore, GtUword l_width,
    gt_rchc_score_t l_minscore)
{
  return (gt_randomcodes_hpcorrect_same_rext_letters(sdata, u_kmerpos, -1,
        t_kmerpos, 1)) &&
    (gt_randomcodes_hpcorrect_compute_hscore_manhattan(NULL,
     GT_RANDOMCODES_HPCORRECT_HMER(sdata, u_sfx) + sdata->k,
     GT_RANDOMCODES_HPCORRECT_HMER(sdata, t_sfx) + sdata->k + 2UL,
     (GtUword)sdata->r - 1UL) >= r_minscore) &&
    (gt_randomcodes_hpcorrect_compute_hscore_manhattan(NULL,
     GT_RANDOMCODES_HPCORRECT_HMER(sdata, u_sfx) + sdata->k - l_width - 1UL,
     GT_RANDOMCODES_HPCORRECT_HMER(sdata, t_sfx) + sdata->k - l_width - 1UL,
     l_width) >= l_minscore);
}

static void gt_randomcodes_hpcorrect_firstpass_correct(
    const GtSeqnumrelpos *snrp, const GtUword *suffixes,
    GtUword nofsuffixes, GtRandomcodesHpcorrectData *sdata,
    GtUword u, GtUword t)
{
  char t_char;
  uint8_t t_hlen, t_hlen_b, t_hlen_a, t_hlen_sum, u_hlen_b;
  GtUword u_sfx, u_relpos, u_seqnum, u_startpos, u_kmerpos, u_kplus1pos,
                t_sfx, t_relpos, t_seqnum, t_startpos, t_kmerpos;
  bool rc;
#ifdef GT_RANDOMCODES_HPCORRECT_DEBUG
  gt_file_xprintf(sdata->outfile,
      "# firstpass_correct, u: "GT_WU", t: "GT_WU"\n", u, t);
#endif
  for (t_sfx = 0; t_sfx < nofsuffixes && sdata->rextset[t_sfx] != t; t_sfx++)
    { /* nothing */ }
  gt_assert(t_sfx < nofsuffixes);
  t_relpos = gt_seqnumrelpos_decode_relpos(snrp, suffixes[t_sfx]);
  t_seqnum = gt_seqnumrelpos_decode_seqnum(snrp, suffixes[t_sfx]);
  t_startpos = gt_encseq_seqstartpos(sdata->encseq, t_seqnum);
  t_kmerpos = t_startpos + t_relpos;
  t_char = gt_encseq_get_decoded_char(sdata->encseq, t_kmerpos + sdata->k,
        GT_READMODE_FORWARD);
  t_hlen_b = GT_RANDOMCODES_HPCORRECT_HMER(sdata, t_sfx)[sdata->k - 1] + 1;
  t_hlen = GT_RANDOMCODES_HPCORRECT_HMER(sdata, t_sfx)[sdata->k] + 1;
  t_hlen_a = GT_RANDOMCODES_HPCORRECT_HMER(sdata, t_sfx)[sdata->k + 1] + 1;
  t_hlen_sum = t_hlen_b + t_hlen + t_hlen_a;
  for (u_sfx = 0; u_sfx < nofsuffixes; u_sfx++)
  {
    if (sdata->rextset[u_sfx] != u)
      continue;
    u_relpos = gt_seqnumrelpos_decode_relpos(snrp, suffixes[u_sfx]);
    u_seqnum = gt_seqnumrelpos_decode_seqnum(snrp, suffixes[u_sfx]);
    u_startpos = gt_encseq_seqstartpos(sdata->encseq, u_seqnum);
    u_kmerpos = u_startpos + u_relpos;
    u_kplus1pos = u_relpos + sdata->k;
    u_hlen_b = GT_RANDOMCODES_HPCORRECT_HMER(sdata, u_sfx)[sdata->k - 1] + 1;
    rc = (u_seqnum >= sdata->nofreads);
    if (rc)
    {
      if (sdata->skip_rc)
        continue;
      u_seqnum = sdata->mirror_nofseqs - 1UL - u_seqnum;
      u_kplus1pos = gt_encseq_seqlength(sdata->encseq, u_seqnum) - 1UL -
        u_kplus1pos;
    }
    if (gt_randomcodes_hpcorrect_similar_rext(sdata, u_sfx, u_kmerpos, 0,
          t_sfx, t_kmerpos, 1, sdata->rext_I_minscore))
    {
      /* deletion found in u */
      sdata->rextset_size[u]--;
      sdata->rextset[u_sfx] = GT_UNDEF_UWORD;
      if (rc)
        gt_file_xprintf(sdata->outfile,""GT_WU"\t"GT_WU"\tI\t%c\t%u\t%u\t%u\n",
            u_seqnum, u_kplus1pos, gt_randomcodes_complement(t_char),
            t_hlen_a, t_hlen, t_hlen_b);
      else
        gt_file_xprintf(sdata->outfile,""GT_WU"\t"GT_WU"\tI\t%c\t%u\t%u\t%u\n",
            u_seqnum, u_kplus1pos - 1UL, t_char,
            t_hlen_b, t_hlen, t_hlen_a);
    }
    else if (gt_randomcodes_hpcorrect_similar_rext(sdata, u_sfx, u_kmerpos, 1,
          t_sfx, t_kmerpos, 1, sdata->rext_R_minscore))
    {
      /* substitution found */
      sdata->rextset_size[u]--;
      sdata->rextset[u_sfx] = GT_UNDEF_UWORD;
      if (rc)
        gt_file_xprintf(sdata->outfile,""GT_WU"\t"GT_WU"\tR\t%c\t0\t%u\t%u\n",
            u_seqnum, u_kplus1pos, gt_randomcodes_complement(t_char),
            t_hlen, t_hlen_b);
      else
        gt_file_xprintf(sdata->outfile,""GT_WU"\t"GT_WU"\tR\t%c\t0\t%u\t%u\n",
            u_seqnum, u_kplus1pos, t_char, t_hlen, t_hlen_a);
    }
    else if (((t_hlen_sum >= u_hlen_b &&
               t_hlen_sum - u_hlen_b <= (uint8_t)1U) ||
              (t_hlen_sum < u_hlen_b &&
               u_hlen_b - t_hlen_sum == (uint8_t)1U)) &&
        gt_randomcodes_hpcorrect_detect_double_insertion(sdata, u_sfx,
          u_kmerpos, t_sfx, t_kmerpos, sdata->rext_J_minscore,
          sdata->rext_J_lwidth, sdata->rext_J_lminscore))
    {
      /* double deletion found in u */
      sdata->rextset_size[u]--;
      sdata->rextset[u_sfx] = GT_UNDEF_UWORD;
      if (rc)
        gt_file_xprintf(sdata->outfile,""GT_WU"\t"GT_WU"\tJ\t%c\t%u\t%u\t%u\n",
            u_seqnum, u_kplus1pos + 1UL, gt_randomcodes_complement(t_char),
            t_hlen_a, t_hlen, t_hlen_b);
      else
        gt_file_xprintf(sdata->outfile,""GT_WU"\t"GT_WU"\tJ\t%c\t%u\t%u\t%u\n",
            u_seqnum, u_kplus1pos - 1UL, t_char,
            t_hlen_b, t_hlen, t_hlen_a);
    }
    else if (gt_randomcodes_hpcorrect_similar_rext(sdata, u_sfx, u_kmerpos, 1,
          t_sfx, t_kmerpos, 0, sdata->rext_D_minscore))
    {
      /* insertion found in u */
      sdata->rextset_size[u]--;
      sdata->rextset[u_sfx] = GT_UNDEF_UWORD;
      gt_file_xprintf(sdata->outfile,""GT_WU"\t"GT_WU"\tD\t-\t0\t0\t%u\n",
                      u_seqnum, u_kplus1pos, t_hlen_a);
    }
  }
}

static inline GtUword
  gt_randomcodes_hpcorrect_firstpass_cluster_right_extensions(
    const GtSeqnumrelpos *snrp, const GtUword *suffixes,
    GtUword nofsuffixes, GtRandomcodesHpcorrectData *sdata,
    GtUword clusternum)
{
  GtUword current_rextset = 0;
  GtUword i, relpos_i, seqnum_i, startpos_i, kmerpos_i,
                j, relpos_j, seqnum_j, startpos_j, kmerpos_j;
#ifdef GT_RANDOMCODES_HPCORRECT_DEBUG
  gt_file_xprintf(sdata->outfile, "# clusternum: "GT_WU", clustersize: "
                  GT_WU "\n", clusternum, sdata->cluster_size[clusternum]);
#endif
  for (i = 0; i < nofsuffixes; i++)
  {
    sdata->rextset[i] = GT_UNDEF_UWORD;
    sdata->rextset_size[i] = 0;
  }
  for (i = 0; i < nofsuffixes; i++)
  {
    if (GT_RANDOMCODES_HPCORRECT_CLUSTERNUM(sdata, i) != clusternum)
      continue;
    if (sdata->rextset[i] != GT_UNDEF_UWORD)
      continue;
    relpos_i = gt_seqnumrelpos_decode_relpos(snrp, suffixes[i]);
    seqnum_i = gt_seqnumrelpos_decode_seqnum(snrp, suffixes[i]);
    startpos_i = gt_encseq_seqstartpos(sdata->encseq, seqnum_i);
    kmerpos_i = startpos_i + relpos_i;
    for (j = i + 1; j < nofsuffixes; j++)
    {
      if (GT_RANDOMCODES_HPCORRECT_CLUSTERNUM(sdata, j) != clusternum)
        continue;
      if (sdata->rextset[j] != GT_UNDEF_UWORD)
        continue;
      relpos_j = gt_seqnumrelpos_decode_relpos(snrp, suffixes[j]);
      seqnum_j = gt_seqnumrelpos_decode_seqnum(snrp, suffixes[j]);
      startpos_j = gt_encseq_seqstartpos(sdata->encseq, seqnum_j);
      kmerpos_j = startpos_j + relpos_j;
      if (gt_randomcodes_hpcorrect_similar_rext(sdata, i, kmerpos_i, 0,
          j, kmerpos_j, 0, sdata->rext_cl_minscore))
      {
        sdata->rextset[j] = current_rextset;
        sdata->rextset_size[current_rextset]++;
      }
    }
    sdata->rextset[i] = current_rextset;
    sdata->rextset_size[current_rextset]++;
    current_rextset++;
  }
  gt_assert(current_rextset > 0);
  return current_rextset;
}

static inline void gt_randomcodes_hpcorrect_correct_hmers(
    const GtSeqnumrelpos *snrp, const GtUword *suffixes,
    GtUword nofsuffixes, GtRandomcodesHpcorrectData *sdata,
    GtUword clusternum)
{
  GtUword i, pos, relpos, seqnum, hpos, maxpos;
  uint8_t *values, consensus, value;
  bool rc = false;
  unsigned int *hdistri, maxuntrust =
     (unsigned int)sdata->cluster_size[clusternum] * sdata->maxuntrustp / 100U;
  for (i = 0; i < nofsuffixes; i++)
  {
    if (GT_RANDOMCODES_HPCORRECT_CLUSTERNUM(sdata, i) != clusternum)
      continue;
    relpos = gt_seqnumrelpos_decode_relpos(snrp, suffixes[i]);
    seqnum = gt_seqnumrelpos_decode_seqnum(snrp, suffixes[i]);
    maxpos = gt_encseq_seqlength(sdata->encseq, seqnum) - 1;
    rc = (seqnum >= sdata->nofreads);
    if (rc)
    {
      if (sdata->skip_rc)
        continue;
      seqnum = sdata->mirror_nofseqs - 1 - seqnum;
    }
    values = GT_RANDOMCODES_HPCORRECT_HMER(sdata, i);
    for (pos = sdata->frompos; pos <= sdata->topos; pos++)
    {
      hpos = rc ? maxpos - (relpos + pos) : relpos + pos;
      if (sdata->skip_read_ends && (hpos == 0 || hpos == maxpos))
        continue;
      consensus = sdata->consensus_hmer[pos];
      if (consensus == GT_HPLSTORE_UNDEF)
        continue;
      value = values[pos];
      if (consensus == value)
        continue;
      hdistri = sdata->hdistri + pos * GT_HPLSTORE_RANGE;
      if (hdistri[value] > maxuntrust)
        continue;
      if (consensus > value)
      {
        gt_file_xprintf(sdata->outfile, ""GT_WU"\t"GT_WU"\tI\t%u\n", seqnum,
                        hpos, (unsigned int)(consensus - value));
      }
      else
      {
        gt_assert(consensus < value);
        gt_file_xprintf(sdata->outfile, ""GT_WU"\t"GT_WU"\tD\t%u\n", seqnum,
                        hpos, (unsigned int)(value - consensus));
      }
    }
  }
}

static inline void gt_randomcode_realloc_hmers_data(GtUword nofsuffixes,
    GtRandomcodesHpcorrectData *sdata)
{
  sdata->hmers_alloc = nofsuffixes + GT_RANDOMCODES_HPCORRECT_HMERS_EXTRA;
  gt_log_log("realloc hmers array to "GT_WU" elements", sdata->hmers_alloc);
  sdata->hmers = gt_realloc(sdata->hmers,
      sizeof (sdata->hmers) * sdata->hmers_alloc * sdata->hmers_width);
  sdata->cluster_size = gt_realloc(sdata->cluster_size,
      sizeof (sdata->cluster_size) * sdata->hmers_alloc);
  if (sdata->union_find == NULL)
  {
    sdata->hmer_cluster = gt_realloc(sdata->hmer_cluster,
        sizeof (sdata->hmer_cluster) * sdata->hmers_alloc);
  }
  else
  {
    sdata->skip = gt_realloc(sdata->skip,
        sizeof (sdata->skip) * sdata->hmers_alloc);
  }
  if (sdata->pw_scores != NULL)
  {
    sdata->pw_scores = gt_realloc(sdata->pw_scores,
        sizeof (sdata->pw_scores) *
          GT_STRICT_TRIANGULAR_SIZE(sdata->hmers_alloc));
    gt_assert(sdata->pw_scores_copy != NULL);
    sdata->pw_scores_copy = gt_realloc(sdata->pw_scores_copy,
        sizeof (sdata->pw_scores_copy) *
          GT_STRICT_TRIANGULAR_SIZE(sdata->hmers_alloc));
  }
  if (sdata->firstpass)
  {
    sdata->rextset = gt_realloc(sdata->rextset,
        sizeof (sdata->rextset) * sdata->hmers_alloc);
    sdata->rextset_size = gt_realloc(sdata->rextset_size,
        sizeof (sdata->rextset_size) * sdata->hmers_alloc);
  }
}

static void gt_randomcodes_hpcorrect_fill_hmers(
    const GtSeqnumrelpos *snrp, const GtUword *suffixes,
    GtUword nofsuffixes, GtRandomcodesHpcorrectData *sdata)
{
  GtUword i;
  for (i = 0; i < nofsuffixes; i++)
  {
    GtUword relpos, seqnum, startpos, kmerpos;
    relpos = gt_seqnumrelpos_decode_relpos(snrp, suffixes[i]);
    seqnum = gt_seqnumrelpos_decode_seqnum(snrp, suffixes[i]);
    startpos = gt_encseq_seqstartpos(sdata->encseq, seqnum);
    kmerpos = startpos + relpos;
    gt_hplstore_get_range(sdata->hplstore,
        GT_RANDOMCODES_HPCORRECT_HMER(sdata,i), kmerpos, sdata->hmers_width);
  }
}

static GtUword gt_randomcodes_hpcorrect_cluster(bool *allidentical,
    GtUword nofsuffixes, GtRandomcodesHpcorrectData *sdata)
{
  GtUword nofclusters;
  gt_assert(allidentical != NULL);
  if (sdata->pw_scores != NULL)
  {
    size_t nofpairs = GT_STRICT_TRIANGULAR_SIZE(nofsuffixes),
           quantile_pos = (nofpairs * sdata->quantile) / 100UL;
    if (quantile_pos > 0)
      quantile_pos -= 1;
    gt_assert(quantile_pos < nofpairs);
    gt_randomcodes_hpcorrect_fill_pw_scores_matrix(nofsuffixes, sdata);
    memcpy(sdata->pw_scores_copy, sdata->pw_scores,
        sizeof (sdata->pw_scores) * nofpairs);
    qsort(sdata->pw_scores_copy, nofpairs, sizeof (sdata->pw_scores_copy),
          gt_rchc_score_t_cmp);
    sdata->clustering_minscore = sdata->pw_scores_copy[quantile_pos];
  }
  nofclusters = sdata->union_find != NULL
    ? gt_randomcodes_hpcorrect_cluster_union_find(allidentical, nofsuffixes,
        sdata->clustering_minscore, sdata)
    : gt_randomcodes_hpcorrect_cluster_greedy(allidentical, nofsuffixes,
        sdata->clustering_minscore, sdata);
#ifdef GT_RANDOMCODES_HPCORRECT_VERBOSE
    gt_file_xprintf(sdata->outfile,
        "# nsfx:\t"GT_WU"\tncls:\t"GT_WU"\n", nofsuffixes, nofclusters);
#endif
#ifdef GT_RANDOMCODES_HPCORRECT_DEBUG
  if (sdata->pw_scores != NULL)
    gt_randomcodes_hpcorrect_show_pw_scores(nofsuffixes, sdata);
#endif
  return nofclusters;
}

static inline void gt_randomcodes_hpcorrect_process_kmer_itv(
    const GtSeqnumrelpos *snrp, const GtUword *suffixes,
    GtUword nofsuffixes, GtRandomcodesHpcorrectData *sdata);

static void gt_randomcodes_hpcorrect_partition_kmer_itv(
    const GtSeqnumrelpos *snrp, const GtUword *suffixes,
    GtUword nofsuffixes, GtRandomcodesHpcorrectData *sdata)
{
    GtUword partnum, partsize, offset, nof_parts;
    gt_assert(sdata->maxwidth > 0);
    nof_parts = 1 + nofsuffixes / sdata->maxwidth;
    partsize = nofsuffixes / nof_parts;
    for (partnum = 0, offset = 0; partnum < nof_parts - 1;
         partnum++, offset+=partsize)
    {
      gt_randomcodes_hpcorrect_process_kmer_itv(snrp, suffixes + offset,
          partsize, sdata);
    }
    if (nofsuffixes - offset > 0)
      gt_randomcodes_hpcorrect_process_kmer_itv(snrp, suffixes + offset,
          nofsuffixes - offset, sdata);
}

static inline void gt_randomcodes_hpcorrect_firstpass_process_cluster(
    const GtSeqnumrelpos *snrp, const GtUword *suffixes,
    GtUword nofsuffixes, GtRandomcodesHpcorrectData *sdata,
    GtUword clusternum)
{
  GtUword nofrextsets, u_rextsetnum, t_rextsetnum, untrusted, trusted;
  nofrextsets = gt_randomcodes_hpcorrect_firstpass_cluster_right_extensions(
      snrp, suffixes, nofsuffixes, sdata, clusternum);
#ifdef GT_RANDOMCODES_HPCORRECT_VERBOSE
  gt_file_xprintf(sdata->outfile, "# nofrextsets:"GT_WU"\n", nofrextsets);
#endif
#ifdef GT_RANDOMCODES_HPCORRECT_DEBUG
  gt_randomcodes_hpcorrect_show_kplus1(snrp, suffixes, nofsuffixes, sdata,
      clusternum);
#endif
  if (nofrextsets == 1UL)
    return;
  untrusted = sdata->cluster_size[clusternum] * sdata->maxuntrustp / 100 + 1UL;
  trusted = sdata->cluster_size[clusternum] * sdata->mintrustp / 100;
  if (trusted <= untrusted)
    trusted = untrusted + 1UL;
#ifdef GT_RANDOMCODES_HPCORRECT_VERBOSE
  gt_file_xprintf(sdata->outfile, "# untrusted: "GT_WU", trusted: "GT_WU"\n",
      untrusted, trusted);
#endif
  for (u_rextsetnum = 0; u_rextsetnum < nofrextsets; u_rextsetnum++)
  {
    if (sdata->rextset_size[u_rextsetnum] <= (unsigned int)untrusted)
    {
#ifdef GT_RANDOMCODES_HPCORRECT_DEBUG
      gt_file_xprintf(sdata->outfile, "# rextset "GT_WU" is untrusted\n",
          u_rextsetnum);
#endif
      for (t_rextsetnum = 0; t_rextsetnum < nofrextsets; t_rextsetnum++)
      {
        if (u_rextsetnum != t_rextsetnum &&
            sdata->rextset_size[t_rextsetnum] >= (unsigned int)trusted)
        {
#ifdef GT_RANDOMCODES_HPCORRECT_DEBUG
          gt_file_xprintf(sdata->outfile, "# ... and rextset " GT_WU
                          " is trusted\n", t_rextsetnum);
#endif
          gt_randomcodes_hpcorrect_firstpass_correct(snrp, suffixes,
              nofsuffixes, sdata, u_rextsetnum, t_rextsetnum);
          /* corrected instances are removed from the rextset */
          if (sdata->rextset_size[u_rextsetnum] == 0)
            break;
        }
      }
    }
  }
}

static inline void gt_randomcodes_hpcorrect_process_cluster(
    const GtSeqnumrelpos *snrp, const GtUword *suffixes,
    GtUword nofsuffixes, GtRandomcodesHpcorrectData *sdata,
    GtUword clusternum)
{
  gt_randomcodes_hpcorrect_compute_consensus_hmer(nofsuffixes, sdata,
      clusternum);
  gt_randomcodes_hpcorrect_correct_hmers(snrp, suffixes, nofsuffixes,
      sdata, clusternum);
}

static inline void gt_randomcodes_hpcorrect_process_kmer_itv(
    const GtSeqnumrelpos *snrp, const GtUword *suffixes,
    GtUword nofsuffixes, GtRandomcodesHpcorrectData *sdata)
{
  GtUword clusternum, nofclusters;
  bool allidentical;
  if (nofsuffixes > sdata->maxwidth)
  {
    gt_randomcodes_hpcorrect_partition_kmer_itv(snrp, suffixes,
        nofsuffixes, sdata);
    return;
  }
#ifdef GT_RANDOMCODES_HPCORRECT_VERBOSE
  gt_file_xprintf(sdata->outfile, "# nofsfx:\t"GT_WU"\n", nofsuffixes);
#endif
  if (nofsuffixes < 3UL)
    return;
  if (sdata->hmers_alloc < nofsuffixes)
    gt_randomcode_realloc_hmers_data(nofsuffixes, sdata);
  gt_randomcodes_hpcorrect_fill_hmers(snrp, suffixes, nofsuffixes, sdata);
  nofclusters = gt_randomcodes_hpcorrect_cluster(&allidentical, nofsuffixes,
      sdata);
  if (allidentical)
    return;
  if (nofclusters > nofsuffixes - 2UL)
    return;
  for (clusternum = 0; clusternum < nofsuffixes; clusternum++)
  {
    if (sdata->cluster_size[clusternum] >= 3UL)
    {
      if (!sdata->firstpass)
        gt_randomcodes_hpcorrect_process_cluster(snrp, suffixes, nofsuffixes,
            sdata, clusternum);
      else
        gt_randomcodes_hpcorrect_firstpass_process_cluster(snrp, suffixes,
            nofsuffixes, sdata, clusternum);
    }
  }
}

int gt_randomcodes_hpcorrect_process_bucket(void *data,
    const GtUword *bucketofsuffixes, const GtSeqnumrelpos *snrp,
    const uint16_t *lcptab_bucket, GtUword numberofsuffixes,
    GT_UNUSED unsigned int sortingdepth, GT_UNUSED GtError *err)
{
  GtUword itvstart, next_itvstart;
  unsigned int lcpvalue;
  GtRandomcodesHpcorrectData *sdata = data;

#ifdef GT_RANDOMCODES_HPCORRECT_DEBUG
  gt_file_xprintf(sdata->outfile, "# bucketsize:\t"GT_WU"\n", numberofsuffixes);
#endif
  for (itvstart = 0, next_itvstart = 1UL; next_itvstart < numberofsuffixes;
      next_itvstart++)
  {
    lcpvalue = (unsigned int) lcptab_bucket[next_itvstart];
#ifdef GT_RANDOMCODES_HPCORRECT_DEBUG
    gt_randomcodes_hpcorrect_show_expanded_kmer(sdata, bucketofsuffixes, snrp,
        next_itvstart - 1UL, lcpvalue);
#endif
    if (lcpvalue < (unsigned int)sdata->k)
    {
#ifdef GT_RANDOMCODES_HPCORRECT_DEBUG
      gt_file_xprintf(sdata->outfile, "# orig-nofsfx:\t"GT_WU"\n",
          next_itvstart - itvstart);
#endif
      gt_randomcodes_hpcorrect_process_kmer_itv(snrp,
          bucketofsuffixes + itvstart, next_itvstart - itvstart, data);
      itvstart = next_itvstart;
    }
  }
#ifdef GT_RANDOMCODES_HPCORRECT_DEBUG
  gt_randomcodes_hpcorrect_show_expanded_kmer(sdata, bucketofsuffixes, snrp,
      next_itvstart - 1UL, lcpvalue);
  gt_file_xprintf(sdata->outfile, "# orig-nofsfx:\t"GT_WU"\n",
      next_itvstart - itvstart);
#endif
  gt_randomcodes_hpcorrect_process_kmer_itv(snrp,
      bucketofsuffixes + itvstart, next_itvstart - itvstart, data);
  return 0;
}

/*
  skip_hmer_ends allows one to avoid wrong corrections like:
     read1    : GAAAaAAAT <- has a mismatch
     consensus: GAAACAAAT
     false corrections: read1 -> GAAAT...
                        read2 -> GAAAaAAACAAAT...
  or:
     read1:     GAAA-AAAT <- an hpol is completely missing
     consensus: GAAACAAAT
     false corrections: same as before

  skip_hmer_ends implies skip_read_ends as read ends are always hmer ends
*/

#define GT_RANDOMCODES_HPCORRECT_FILESFX ".cor"
#define GT_RANDOMCODES_HPCORRECT_FIRSTPASS_FILESFX ".fpc"
#define GT_RANDOMCODES_HPCORRECT_UNION_FIND_INIT 1024UL

GtRandomcodesHpcorrectData *gt_randomcodes_hpcorrect_data_new(
    GtEncseq *encseq, GtHplstore *hplstore, unsigned int k,
    bool firstpass, unsigned int r, unsigned int mintrustp,
    unsigned int maxuntrustp, bool greedy_clustering,
    bool skip_read_ends, bool skip_hmer_ends, bool skip_rc,
    bool non_redundant, bool best_score_clustering, bool manhattan,
    GtWord clustering_param, GtUword maxwidth, int rext_cl_minscore,
    int rext_I_minscore, int rext_J_minscore, int rext_R_minscore,
    int rext_D_minscore, int rext_J_lminscore, GtUword rext_J_lwidth,
    GtStr *indexname, unsigned int threadnum, GtError *err)
{
  GtRandomcodesHpcorrectData *sdata = gt_malloc(sizeof *sdata);
  GtStr *path = gt_str_clone(indexname);
  gt_str_append_char(path, '.');
  gt_str_append_uint(path, threadnum);
  if (firstpass)
    gt_str_append_cstr(path, GT_RANDOMCODES_HPCORRECT_FIRSTPASS_FILESFX);
  else
    gt_str_append_cstr(path, GT_RANDOMCODES_HPCORRECT_FILESFX);
  sdata->k = (GtUword)k;
  sdata->firstpass = firstpass;
  sdata->mintrustp = mintrustp;
  sdata->maxuntrustp = maxuntrustp;
  sdata->maxwidth = maxwidth;
  sdata->manhattan = manhattan;
  sdata->skip_read_ends = skip_read_ends;
  sdata->rext_cl_minscore = (gt_rchc_score_t)rext_cl_minscore;
  sdata->rext_I_minscore = (gt_rchc_score_t)rext_I_minscore;
  sdata->rext_J_minscore = (gt_rchc_score_t)rext_J_minscore;
  sdata->rext_R_minscore = (gt_rchc_score_t)rext_R_minscore;
  sdata->rext_D_minscore = (gt_rchc_score_t)rext_D_minscore;
  sdata->rext_J_lminscore = (gt_rchc_score_t)rext_J_lminscore;
  sdata->rext_J_lwidth = rext_J_lwidth;
  if (skip_hmer_ends)
  {
    sdata->frompos = 1UL;
    gt_assert(sdata->k > 2UL);
    sdata->topos = sdata->k - 2UL;
  }
  else
  {
    sdata->frompos = 0;
    sdata->topos = sdata->k - 1UL;
  }
  if (non_redundant)
    sdata->frompos = sdata->topos;
  sdata->skip_rc = skip_rc;
  sdata->outfile = gt_file_new(gt_str_get(path), "w", err);
  gt_str_delete(path);
  sdata->encseq = encseq;
  sdata->hplstore = hplstore;
  sdata->mirror_nofseqs = gt_encseq_num_of_sequences(encseq);
  sdata->nofreads = sdata->mirror_nofseqs;
  if (gt_encseq_is_mirrored(sdata->encseq))
    sdata->nofreads >>= 1;
  sdata->hdistri = gt_malloc(sizeof (sdata->hdistri) * GT_HPLSTORE_RANGE *
      sdata->k);
  sdata->hmers_alloc = GT_RANDOMCODES_HPCORRECT_HMERS_INIT;
  if (best_score_clustering)
  {
    gt_assert(clustering_param >= 0 && clustering_param <= 100L);
    sdata->quantile = (unsigned int)clustering_param;
    sdata->pw_scores = gt_malloc(sizeof (sdata->pw_scores) *
            GT_STRICT_TRIANGULAR_SIZE(sdata->hmers_alloc));
    sdata->pw_scores_copy = gt_malloc(sizeof (sdata->pw_scores_copy) *
            GT_STRICT_TRIANGULAR_SIZE(sdata->hmers_alloc));
  }
  else
  {
    sdata->clustering_minscore = clustering_param;
    sdata->pw_scores = NULL;
    sdata->pw_scores_copy = NULL;
  }
  sdata->cluster_size = gt_malloc(sizeof (sdata->cluster_size) *\
      sdata->hmers_alloc);
  sdata->consensus_hmer = gt_malloc(sizeof (sdata->consensus_hmer) * sdata->k);
  if (greedy_clustering)
  {
    sdata->hmer_cluster = gt_malloc(sizeof (sdata->hmer_cluster) *
        sdata->hmers_alloc);
    sdata->union_find = NULL;
    sdata->skip = NULL;
  }
  else
  {
    sdata->union_find =
      gt_union_find_new(GT_RANDOMCODES_HPCORRECT_UNION_FIND_INIT);
    sdata->hmer_cluster = NULL;
    sdata->skip = gt_malloc(sizeof (sdata->skip) * sdata->hmers_alloc);
  }
  if (sdata->firstpass)
  {
    sdata->r = r;
    sdata->rextset = gt_malloc(sizeof (sdata->rextset) * sdata->hmers_alloc);
    sdata->rextset_size = gt_malloc(sizeof (sdata->rextset_size) *
        sdata->hmers_alloc);
    sdata->esr1 = gt_encseq_create_reader_with_readmode(sdata->encseq,
        GT_READMODE_FORWARD, 0);
    sdata->esr2 = gt_encseq_create_reader_with_readmode(sdata->encseq,
        GT_READMODE_FORWARD, 0);
    sdata->hmers_width = sdata->k + sdata->r + 1UL;
  }
  else
  {
    sdata->rextset = NULL;
    sdata->rextset_size = NULL;
    sdata->esr1 = NULL;
    sdata->esr2 = NULL;
    sdata->hmers_width = sdata->k;
  }
  sdata->hmers = gt_malloc(sizeof (sdata->hmers) * sdata->hmers_alloc *
      sdata->hmers_width);
  return sdata;
}

void gt_randomcodes_hpcorrect_data_delete(GtRandomcodesHpcorrectData *sdata)
{
  gt_assert(sdata != NULL);
  gt_free(sdata->hmers);
  gt_free(sdata->pw_scores);
  gt_free(sdata->pw_scores_copy);
  gt_free(sdata->hmer_cluster);
  gt_free(sdata->cluster_size);
  gt_free(sdata->skip);
  gt_free(sdata->consensus_hmer);
  gt_free(sdata->hdistri);
  gt_union_find_delete(sdata->union_find);
  gt_free(sdata->rextset);
  gt_free(sdata->rextset_size);
  gt_encseq_reader_delete(sdata->esr1);
  gt_encseq_reader_delete(sdata->esr2);
  gt_file_delete(sdata->outfile);
  gt_free(sdata);
}
