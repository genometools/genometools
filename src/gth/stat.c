/*
  Copyright (c) 2004-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2004-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/ma_api.h"
#include "core/safearith.h"
#include "gth/default.h"
#include "gth/gthoutput.h"
#include "gth/time.h"
#include "gth/stat.h"

#define INFOCHAR\
        if (!stat->gthfilestat_mode) {          \
          gt_file_xprintf(outfp, "%c ", COMMENTCHAR); \
        }

struct GthStat {
  bool exondistri,       /* compute/output exon length distribution */
       introndistri,     /* compute/output intron length distribution */
       matchnumdistri,   /* compute/output match number distribution */
       refseqcovdistri,  /* compute/output reference sequence coverage
                            distribution */
       sa_stats,         /* compute/output spliced alignment statistics (e.g.,
                            alignment score and coverage distributions) */
       gthfilestat_mode; /* modify output for gthfilestat */

  unsigned long numofchains,             /* total number of computed chains */
       numofremovedzerobaseexons,        /* see function replacezerobaseexons()
                                            for explanation */
       numofautointroncutoutcalls,       /* number of times when the intron
                                            cutout technique was used
                                            automatically */
       numofunsuccessfulintroncutoutDPs, /* corresponds to
                                            ERROR_CUTOUT_NOT_IN_INTRON */
       numoffailedDPparameterallocations,/* correspoinds to
                                            ERROR_DP_PARAMETER_ALLOCATION_FAILED
                                          */
       numoffailedmatrixallocations,     /* corresponds to
                                            ERROR_MATRIX_ALLOCATION_FAILED */
       numofundeterminedSAs,             /* corresponds to
                                            ERROR_SA_COULD_NOT_BE_DETERMINED */
       numoffilteredpolyAtailmatches,    /* the number of filtered out matches
                                            which represent a matching poly(A)
                                            tail */

  /* memory statistics */
       numofSAs,                         /* number of computed SAs */
       numofPGLs_stored,                 /* number of stored PGLs */
       totalsizeofbacktracematricesinMB,
       numofbacktracematrixallocations;

  /* distributions */
  GtDiscDistri *exondistribution,
               *introndistribution,
               *matchnumdistribution,
               *refseqcoveragedistribution,
               *sa_alignment_score_distribution,
               *sa_coverage_distribution;
};

GthStat *gth_stat_new(void)
{
  GthStat *stat;

  stat = gt_malloc(sizeof (GthStat));

  stat->exondistri       = GTH_DEFAULT_EXONDISTRI;
  stat->introndistri     = GTH_DEFAULT_INTRONDISTRI;
  stat->matchnumdistri   = GTH_DEFAULT_MATCHNUMDISTRI;
  stat->refseqcovdistri  = GTH_DEFAULT_REFSEQCOVDISTRI;
  stat->sa_stats         = false;
  stat->gthfilestat_mode = false;

  stat->numofchains                       = 0;
  stat->numofremovedzerobaseexons         = 0;
  stat->numofautointroncutoutcalls        = 0;
  stat->numofunsuccessfulintroncutoutDPs  = 0;
  stat->numoffailedDPparameterallocations = 0;
  stat->numoffailedmatrixallocations      = 0;
  stat->numofundeterminedSAs              = 0;
  stat->numoffilteredpolyAtailmatches     = 0;

  /* init variables for memory statistics */
  stat->numofSAs                          = 0;
  stat->numofPGLs_stored                  = 0;
  stat->totalsizeofbacktracematricesinMB  = 0;
  stat->numofbacktracematrixallocations   = 0;

  /* init distributions */
  stat->exondistribution = gt_disc_distri_new();
  stat->introndistribution = gt_disc_distri_new();
  stat->matchnumdistribution = gt_disc_distri_new();
  stat->refseqcoveragedistribution = gt_disc_distri_new();
  stat->sa_alignment_score_distribution = gt_disc_distri_new();
  stat->sa_coverage_distribution = gt_disc_distri_new();

  return stat;
}

void gth_stat_enable_exondistri(GthStat *stat)
{
  gt_assert(stat);
  stat->exondistri = true;
}

void gth_stat_enable_introndistri(GthStat *stat)
{
  gt_assert(stat);
  stat->introndistri= true;
}

void gth_stat_enable_matchnumdistri(GthStat *stat)
{
  gt_assert(stat);
  stat->matchnumdistri = true;
}

void gth_stat_enable_refseqcovdistri(GthStat *stat)
{
  gt_assert(stat);
  stat->refseqcovdistri = true;
}

void gth_stat_enable_sa_stats(GthStat *stat)
{
  gt_assert(stat);
  stat->sa_stats = true;
}

void gth_stat_enable_gthfilestat_mode(GthStat *stat)
{
  gt_assert(stat);
  stat->gthfilestat_mode= true;
}

void gth_stat_increment_numofunsuccessfulintroncutoutDPs(GthStat *stat)
{
  gt_assert(stat);
  stat->numofunsuccessfulintroncutoutDPs++;
}

void gth_stat_increment_numofundeterminedSAs(GthStat *stat)
{
  gt_assert(stat);
  stat->numofundeterminedSAs++;
}

void gth_stat_increment_numofautointroncutoutcalls(GthStat *stat)
{
  gt_assert(stat);
  stat->numofautointroncutoutcalls++;
}

void gth_stat_increment_numoffailedmatrixallocations(GthStat *stat)
{
  gt_assert(stat);
  stat->numoffailedmatrixallocations++;
}

void gth_stat_increment_numoffailedDPparameterallocations(GthStat *stat)
{
  gt_assert(stat);
  stat->numoffailedDPparameterallocations++;
}

void gth_stat_increment_numofbacktracematrixallocations(GthStat *stat)
{
  gt_assert(stat);
  stat->numofbacktracematrixallocations++;
}

void gth_stat_increment_numofremovedzerobaseexons(GthStat *stat)
{
  gt_assert(stat);
  stat->numofremovedzerobaseexons++;
}

void gth_stat_increment_numofSAs(GthStat *stat)
{
  gt_assert(stat);
  stat->numofSAs++;
}

void gth_stat_increase_numofchains(GthStat *stat, unsigned long addend)
{
  gt_assert(stat);
  stat->numofchains += addend;
}

void gth_stat_increase_totalsizeofbacktracematricesinMB(GthStat *stat,
                                                        unsigned long addend)
{
  gt_assert(stat);
  gt_safe_add(stat->totalsizeofbacktracematricesinMB,
              stat->totalsizeofbacktracematricesinMB, addend);
}

void gth_stat_increase_numofPGLs_stored(GthStat *stat, unsigned long addend)
{
  gt_assert(stat);
  stat->numofPGLs_stored += addend;
}

unsigned long gth_stat_get_numofSAs(GthStat *stat)
{
  gt_assert(stat);
  return stat->numofSAs;
}

bool gth_stat_get_exondistri(GthStat *stat)
{
  gt_assert(stat);
  return stat->exondistri;
}

bool gth_stat_get_introndistri(GthStat *stat)
{
  gt_assert(stat);
  return stat->introndistri;
}

GtDiscDistri* gth_stat_get_exondistribution(GthStat *stat)
{
  gt_assert(stat);
  return stat->exondistribution;
}

GtDiscDistri* gth_stat_get_introndistribution(GthStat *stat)
{
  gt_assert(stat);
  return stat->introndistribution;
}

bool gth_stat_get_matchnumdistri(GthStat *stat)
{
  gt_assert(stat);
  return stat->matchnumdistri;
}

bool gth_stat_get_refseqcovdistri(GthStat *stat)
{
  gt_assert(stat);
  return stat->refseqcovdistri;
}

void gth_stat_add_to_matchnumdistri(GthStat *stat, unsigned long data)
{
  gt_assert(stat);
  if (stat->matchnumdistri)
    gt_disc_distri_add(stat->matchnumdistribution, data);
}

void gth_stat_add_to_refseqcovdistri(GthStat *stat, unsigned long data)
{
  gt_assert(stat);
  if (stat->refseqcovdistri)
    gt_disc_distri_add(stat->refseqcoveragedistribution, data);
}

void gth_stat_add_to_sa_alignment_score_distri(GthStat *stat,
                                               unsigned long data)
{
  gt_assert(stat);
  if (stat->sa_stats)
    gt_disc_distri_add(stat->sa_alignment_score_distribution, data);
}

void gth_stat_add_to_sa_coverage_distri(GthStat *stat, unsigned long data)
{
  gt_assert(stat);
  if (stat->sa_stats)
    gt_disc_distri_add(stat->sa_coverage_distribution, data);
}

static void outputgeneralstatistics(GthStat *stat, bool show_full_stats,
                                    GtFile *outfp)
{
  gt_assert(stat);
  if (show_full_stats) {
    gt_file_xprintf(outfp, "%c general statistics:\n", COMMENTCHAR);
    switch (stat->numofchains) {
      case 0:
        gt_file_xprintf(outfp, "%c no chain has been computed\n",
                           COMMENTCHAR);
        break;
      case 1:
        gt_file_xprintf(outfp, "%c 1 chain has been computed\n",
                           COMMENTCHAR);
        break;
      default:
        gt_file_xprintf(outfp, "%c %lu chains have been computed\n",
                        COMMENTCHAR, stat->numofchains);
    }
  }
}

static void outputmemorystatistics(GthStat *stat, bool show_full_stats,
                                   GtFile *outfp)
{
  gt_assert(stat);

  /* SA memory statistics */
  INFOCHAR;
  gt_file_xprintf(outfp, "memory statistics:\n");
  INFOCHAR;
  gt_file_xprintf(outfp, "%lu spliced alignments have been stored\n",
            stat->numofSAs);

  /* PGL memory statistics */
  INFOCHAR;
  gt_file_xprintf(outfp, "%lu predicted gene locations have been stored\n",
            stat->numofPGLs_stored);

  /* DP memory statistics */
  if (show_full_stats) {
    gt_file_xprintf(outfp,
              "%c %lu megabytes was the average size of the backtrace matrix\n",
              COMMENTCHAR, (stat->numofbacktracematrixallocations > 0
                           ? stat->totalsizeofbacktracematricesinMB /
                             stat->numofbacktracematrixallocations
                           : 0));
    switch (stat->numofbacktracematrixallocations) {
      case 0:
        gt_file_xprintf(outfp, "%c no backtrace matrix has been allocated\n",
                  COMMENTCHAR);
        break;
      case 1:
        gt_file_xprintf(outfp, "%c 1 backtrace matrix has been allocated\n",
                  COMMENTCHAR);
        break;
      default:
        gt_file_xprintf(outfp,"%c %lu backtrace matrices have been "
                        "allocated\n", COMMENTCHAR,
                        stat->numofbacktracematrixallocations);
    }
  }
}

void gth_stat_show(GthStat *stat, bool show_full_stats, bool xmlout,
                   GtFile *outfp)
{
  char *timestring;

  gt_assert(stat);

  /* begin XML comment */
  if (xmlout)
    gt_file_xprintf(outfp, "<!--\n");

  /* output exon length distribution */
  if (stat->exondistri) {
    gt_file_xprintf(outfp, "%c length distribution of all exons:\n",
                    COMMENTCHAR);
    gt_disc_distri_show(stat->exondistribution, outfp);
  }

  /* output intron length distribution */
  if (stat->introndistri) {
    if (stat->exondistri)
      gt_file_xprintf(outfp, "%c\n", COMMENTCHAR);
    gt_file_xprintf(outfp, "%c length distribution of all introns:\n",
                    COMMENTCHAR);
    gt_disc_distri_show(stat->introndistribution, outfp);
  }

  /* output match number distribution */
  if (stat->matchnumdistri) {
    if (stat->exondistri || stat->introndistri)
      gt_file_xprintf(outfp, "%c\n", COMMENTCHAR);
    gt_file_xprintf(outfp, "%c distribution of match numbers (per genomic "
                    "file, per reference sequence:\n", COMMENTCHAR);
    gt_disc_distri_show(stat->matchnumdistribution, outfp);
  }

  /* output reference sequence coverage distribution */
  if (stat->refseqcovdistri) {
    if (stat->exondistri || stat->introndistri || stat->matchnumdistri)
      gt_file_xprintf(outfp, "%c\n", COMMENTCHAR);
    gt_file_xprintf(outfp, "%c reference sequence coverage distribution (of "
                    "global chains):\n", COMMENTCHAR);
    gt_disc_distri_show(stat->refseqcoveragedistribution, outfp);
  }

  /* output spliced alignment statistics */
  if (stat->sa_stats) {
    if (stat->exondistri     || stat->introndistri ||
        stat->matchnumdistri || stat->refseqcovdistri) {
      gt_file_xprintf(outfp, "%c\n", COMMENTCHAR);
    }
    INFOCHAR;
    gt_file_xprintf(outfp,
                       "spliced alignment alignment score distribution:\n");
    gt_disc_distri_show(stat->sa_alignment_score_distribution, outfp);
    INFOCHAR;
    gt_file_xfputc('\n', outfp);
    INFOCHAR;
    gt_file_xprintf(outfp, "spliced alignment coverage distribution:\n");
    gt_disc_distri_show(stat->sa_coverage_distribution, outfp);
  }

  /* output general statistics */
  outputgeneralstatistics(stat, show_full_stats, outfp);
  INFOCHAR;
  gt_file_xfputc('\n', outfp);

  /* output the memory statistics */
  outputmemorystatistics(stat, show_full_stats, outfp);

  /* output time */
  INFOCHAR;
  gt_file_xfputc('\n', outfp);
  INFOCHAR;
  timestring = gth_get_time();
  gt_file_xprintf(outfp, "date finished: %s\n", timestring);
  gt_free(timestring);

  /* output important messages */
  if (stat->numofremovedzerobaseexons         ||
      stat->numofautointroncutoutcalls        ||
      stat->numofunsuccessfulintroncutoutDPs  ||
      stat->numoffailedDPparameterallocations ||
      stat->numoffailedmatrixallocations      ||
      stat->numofundeterminedSAs              ||
      stat->numoffilteredpolyAtailmatches) {
    gt_file_xprintf(outfp, "%c\n", COMMENTCHAR);
    gt_file_xprintf(outfp, "%c important messages:\n", COMMENTCHAR);
    if (stat->numofremovedzerobaseexons > 0) {
      gt_file_xprintf(outfp, "%c %lu removed zero base exons\n",
                         COMMENTCHAR, stat->numofremovedzerobaseexons);
    }
    if (stat->numofautointroncutoutcalls > 0) {
      gt_file_xprintf(outfp, "%c %lu times the intron cutout technique was "
                         "used automatically\n", COMMENTCHAR,
                         stat->numofautointroncutoutcalls);
    }
    if (stat->numofunsuccessfulintroncutoutDPs > 0) {
      gt_file_xprintf(outfp, "%c %lu unsuccessful DP calls using intron "
                         "cutout technique\n", COMMENTCHAR,
                         stat->numofunsuccessfulintroncutoutDPs);
    }
    if (stat->numoffailedDPparameterallocations > 0) {
      gt_file_xprintf(outfp, "%c %lu DP parameter allocations failed\n",
                         COMMENTCHAR, stat->numoffailedDPparameterallocations);
    }
    if (stat->numoffailedmatrixallocations > 0) {
      gt_file_xprintf(outfp, "%c %lu matrix allocations failed\n",
                         COMMENTCHAR, stat->numoffailedmatrixallocations);
    }
    if (stat->numofundeterminedSAs > 0) {
      gt_file_xprintf(outfp, "%c %lu undetermined spliced alignments\n",
                         COMMENTCHAR, stat->numofundeterminedSAs);
    }
    if (stat->numoffilteredpolyAtailmatches > 0) {
      gt_file_xprintf(outfp,
                      "%c %lu matches containing a poly(A) tail filtered\n",
                         COMMENTCHAR, stat->numoffilteredpolyAtailmatches);
    }
  }

  /* end XML comment */
  if (xmlout)
    gt_file_xprintf(outfp, "-->\n");
}

void gth_stat_delete(GthStat *stat)
{
  if (!stat) return;
  gt_disc_distri_delete(stat->exondistribution);
  gt_disc_distri_delete(stat->introndistribution);
  gt_disc_distri_delete(stat->matchnumdistribution);
  gt_disc_distri_delete(stat->refseqcoveragedistribution);
  gt_disc_distri_delete(stat->sa_alignment_score_distribution);
  gt_disc_distri_delete(stat->sa_coverage_distribution);
  gt_free(stat);
}
