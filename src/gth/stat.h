/*
  Copyright (c) 2004-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef STAT_H
#define STAT_H

#include <stdbool.h>
#include "core/disc_distri_api.h"
#include "core/error.h"

/* stores all kinds of gth related statistical information */
typedef struct GthStat GthStat;

GthStat*      gth_stat_new(void);
void          gth_stat_enable_exondistri(GthStat*);
void          gth_stat_enable_introndistri(GthStat*);
void          gth_stat_enable_matchnumdistri(GthStat*);
void          gth_stat_enable_refseqcovdistri(GthStat*);
void          gth_stat_enable_sa_stats(GthStat*);
void          gth_stat_enable_gthfilestat_mode(GthStat*);
void          gth_stat_increment_numofunsuccessfulintroncutoutDPs(GthStat*);
void          gth_stat_increment_numofundeterminedSAs(GthStat*);
void          gth_stat_increment_numofautointroncutoutcalls(GthStat*);
void          gth_stat_increment_numoffailedmatrixallocations(GthStat*);
void          gth_stat_increment_numoffailedDPparameterallocations(GthStat*);
void          gth_stat_increment_numofbacktracematrixallocations(GthStat*);
void          gth_stat_increment_numofremovedzerobaseexons(GthStat*);
void          gth_stat_increment_numofSAs(GthStat*);
void          gth_stat_increase_numofchains(GthStat*, unsigned long);
void          gth_stat_increase_totalsizeofbacktracematricesinMB(GthStat*,
                                                                 unsigned long);
void          gth_stat_increase_numofPGLs_stored(GthStat*, unsigned long);
unsigned long gth_stat_get_numofSAs(GthStat*);
bool          gth_stat_get_exondistri(GthStat*);
bool          gth_stat_get_introndistri(GthStat*);
GtDiscDistri* gth_stat_get_exondistribution(GthStat*);
GtDiscDistri* gth_stat_get_introndistribution(GthStat*);
bool          gth_stat_get_matchnumdistri(GthStat*);
bool          gth_stat_get_refseqcovdistri(GthStat*);
void          gth_stat_add_to_matchnumdistri(GthStat*, unsigned long);
void          gth_stat_add_to_refseqcovdistri(GthStat*, unsigned long);
void          gth_stat_add_to_sa_alignment_score_distri(GthStat*,
                                                        unsigned long);
void          gth_stat_add_to_sa_coverage_distri(GthStat*, unsigned long);
void          gth_stat_show(GthStat*, bool show_full_stats, bool xmlout,
                            GtFile*);
void          gth_stat_delete(GthStat*);

#endif
