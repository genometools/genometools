/*
  Copyright (c) 2003-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef BACKTRACE_PATH_H
#define BACKTRACE_PATH_H

#include "core/file.h"
#include "gth/cutoffmode.h"
#include "gth/gthalignment.h"
#include "gth/gthalphatype.h"
#include "gth/stat.h"

typedef struct GthBacktracePath GthBacktracePath;

#include "gth/path_walker.h"

typedef struct {
  unsigned long genomiccutoff,   /* the number of bp, which will be cut
                                    off when showing the alignment */
                referencecutoff, /* the number of reference bp, which will be
                                    cut off when showing the alignment */
                eopcutoff;       /* the number of Editoperations which will be
                                    cut off when showing the alignment;
                                    this ususally differs from <enomiccutoff>
                                    and <eferencecutoff> because whole introns
                                    are encoded as one multi edit operation */
} Cutoffs;

GthBacktracePath* gth_backtrace_path_new(unsigned long gen_dp_start,
                                         unsigned long gen_dp_length,
                                         unsigned long ref_dp_start,
                                         unsigned long ref_dp_length);
unsigned long   gth_backtrace_path_gen_dp_start(const GthBacktracePath*);
void            gth_backtrace_path_set_gen_dp_start(GthBacktracePath*,
                                                    unsigned long);
unsigned long   gth_backtrace_path_gen_dp_length(const GthBacktracePath*);
void            gth_backtrace_path_set_gen_dp_length(GthBacktracePath*,
                                                     unsigned long);
unsigned long   gth_backtrace_path_ref_dp_length(const GthBacktracePath*);
void            gth_backtrace_path_set_ref_dp_length(GthBacktracePath*,
                                                     unsigned long);
unsigned long   gth_backtrace_path_indelcount(const GthBacktracePath*);
unsigned long   gth_backtrace_path_genomiccutoff_start(const GthBacktracePath*);
unsigned long   gth_backtrace_path_referencecutoff_start(const
                                                         GthBacktracePath*);
unsigned long   gth_backtrace_path_eopcutoff_start(const GthBacktracePath*);
unsigned long   gth_backtrace_path_genomiccutoff_end(const GthBacktracePath*);
unsigned long   gth_backtrace_path_referencecutoff_end(const GthBacktracePath*);
unsigned long   gth_backtrace_path_eopcutoff_end(const GthBacktracePath*);
void            gth_backtrace_path_set_cutoffs_start(GthBacktracePath*,
                                                     Cutoffs*);
void            gth_backtrace_path_set_cutoffs_end(GthBacktracePath*,
                                                   Cutoffs*);
GthAlphatype    gth_backtrace_path_alphatype(const GthBacktracePath*);
void            gth_backtrace_path_set_alphatype(GthBacktracePath*,
                                                 GthAlphatype);
void            gth_backtrace_path_determine_cutoffs(GthBacktracePath*,
                                                     GthCutoffmode
                                                     leadcutoffsmode,
                                                     GthCutoffmode
                                                     termcutoffsmode,
                                                     unsigned long
                                                     cutoffsminexonlen);
void            gth_backtrace_path_remove_zero_base_exons(GthBacktracePath*,
                                                          GthStat*);
bool            gth_backtrace_path_contains_no_zero_base_exons(const
                                                             GthBacktracePath*);

/* add <length> times edit operations of type <eoptype> */
void            gth_backtrace_path_add_eop(GthBacktracePath*, Eoptype eoptype,
                                       unsigned long length);

/* add the corresponding edit operation */
void            gth_backtrace_path_add_match(GthBacktracePath*,
                                             bool ensure_single_match);
void            gth_backtrace_path_add_mismatch(GthBacktracePath*);
void            gth_backtrace_path_add_deletion(GthBacktracePath*);
void            gth_backtrace_path_add_insertion(GthBacktracePath*);
void            gth_backtrace_path_add_intron(GthBacktracePath*);

/* add the corresponding edit operation (only for protein eops!) */
void            gth_backtrace_path_add_mismatch_with_1_gap(GthBacktracePath*);
void            gth_backtrace_path_add_mismatch_with_2_gaps(GthBacktracePath*);
void            gth_backtrace_path_add_deletion_with_1_gap(GthBacktracePath*);
void            gth_backtrace_path_add_deletion_with_2_gaps(GthBacktracePath*);
void            gth_backtrace_path_add_intron_with_1_base_left(
                                                             GthBacktracePath*);
void            gth_backtrace_path_add_intron_with_2_bases_left(
                                                             GthBacktracePath*);
/* adds a dummy which represents a match or a mismatch */
void            gth_backtrace_path_add_dummy(GthBacktracePath*);
/* sets the dummy to a match, if <match> is true; to a mismatch otherwise */
void            gth_backtrace_path_set_dummy(GthBacktracePath*, bool match);
bool            gth_backtrace_path_contain_dummy(const GthBacktracePath*);

bool            gth_backtrace_path_last_is_intron(const GthBacktracePath*);

/* Reverses the edit operations. Necessary, if the editoperations have been
   added in forward direction (reverse direction is the common case) */
void            gth_backtrace_path_reverse(GthBacktracePath*);

/* The following function ensures that the editoperations before introns
   which are not in phase has the maximum length 1.
   Note that the edit operations have to be reversed already, so the
   editoperations checked are actually _after_ the introns not in phase.
   This function makes only sense for protein edit operations, but can be called
   for dna edit operations without any harm */
void            gth_backtrace_path_ensure_length_1_before_introns(
                                                             GthBacktracePath*);

size_t          gth_backtrace_path_sizeof (const GthBacktracePath*);
void            gth_backtrace_path_show(const GthBacktracePath*, bool xmlout,
                                        unsigned int indentlevel, GtFile*);
void            gth_backtrace_path_show_complete(const GthBacktracePath*,
                                                 bool xmlout,
                                                 unsigned int indentlevel,
                                                 GtFile*);
void            gth_backtrace_path_cutoff_start(GthBacktracePath*);
void            gth_backtrace_path_cutoff_end(GthBacktracePath*);
void            gth_backtrace_path_cutoff_walked_path(GthBacktracePath*,
                                                      const GthPathWalker*,
                                                      bool showeops,
                                                      GtFile*);
void            gth_backtrace_path_prepend(GthBacktracePath*,
                                           const GthBacktracePath*);
void            gth_backtrace_path_append(GthBacktracePath*,
                                          const GthBacktracePath*);
bool            gth_backtrace_path_is_valid(const GthBacktracePath*);
void            gth_backtrace_path_reset(GthBacktracePath*);
void            gth_backtrace_path_delete(GthBacktracePath*);

/* XXX: the following functions interface with legacy code, they should be
   removed */

/* return the underlying edit operations (without cutoffs) */
Editoperation*  gth_backtrace_path_get(const GthBacktracePath*);
/* return the number of underlying edit operations (without cutoffs) */
unsigned long   gth_backtrace_path_length(const GthBacktracePath*);

#endif
