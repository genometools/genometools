/*
  Copyright (c) 2014 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2014 Center for Bioinformatics, University of Hamburg

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

#ifndef BLAST_PROCESS_CALL_H
#define BLAST_PROCESS_CALL_H

#include <stdio.h>

#include "core/error_api.h"

/* The <GtBlastProcessCall> class, a helper class to create blast calls. */
typedef struct GtBlastProcessCall GtBlastProcessCall;

/* Returns new <GtBlastProcessCall> object calling blastall for nucleotide
   sequences.
   blastall has to be installed in PATH or GT_BLAST_PATH has to be set to the
   correct path. */
GtBlastProcessCall* gt_blast_process_call_new_all_nucl(void);

/* Returns new <GtBlastProcessCall> object calling blastn.
   blastn has to be installed in PATH or GT_BLAST_PATH has to be set to the
   correct path. */
GtBlastProcessCall* gt_blast_process_call_new_nucl(void);

/* Returns new <GtBlastProcessCall> object calling blastall for protein
   sequences.
   blastall has to be installed in PATH or GT_BLAST_PATH has to be set to the
   correct path. */
GtBlastProcessCall* gt_blast_process_call_new_all_prot(void);

/* Returns new <GtBlastProcessCall> object calling blastp.
   blastp has to be installed in PATH or GT_BLAST_PATH has to be set to the
   correct path. */
GtBlastProcessCall* gt_blast_process_call_new_prot(void);

/* Sets query option for <call> to <query>. See blast documentation for
   explanation of option! */
void                gt_blast_process_call_set_query(GtBlastProcessCall *call,
                                                    const char *query);
/* Sets db option for <call> to <db>. See blast documentation for explanation of
   option! */
void                gt_blast_process_call_set_db(GtBlastProcessCall *call,
                                                 const char *db);
/* Sets evalue option for <call> to <evalue>. See blast documentation for
   explanation of option! */
void                gt_blast_process_call_set_evalue(GtBlastProcessCall *call,
                                                     double evalue);
/* Sets wordsize option for <call> to <wordsize>. See blast documentation for
   explanation of option! */
void                gt_blast_process_call_set_wordsize(GtBlastProcessCall *call,
                                                       int word_size);
/* Sets gapopen option for <call> to <gapopen>. See blast documentation for
   explanation of option! */
void                gt_blast_process_call_set_gapopen(GtBlastProcessCall *call,
                                                      int gapopen);
/* Sets gapextend option for <call> to <gapextend>. See blast documentation for
explanation of option!*/
void                gt_blast_process_call_set_gapextend(
                                                       GtBlastProcessCall *call,
                                                       int gapextend);
/* Sets penalty option for <call> to <penalty>. Nucleotide only option! See
   blast documentation for explanation of option! */
void                gt_blast_process_call_set_penalty(GtBlastProcessCall *call,
                                                      int penalty);
/* Sets reward option for <call> to <reward>. Nucleotide only option! See blast
   documentation for explanation of option! */
void                gt_blast_process_call_set_reward(GtBlastProcessCall *call,
                                                     int reward);
/* Sets number of threads option for <call> to <num_threads>. See blast
   documentation for explanation of option! */
void                gt_blast_process_call_set_num_threads(
                                                       GtBlastProcessCall *call,
                                                       int num_threads);
/* Sets xdrop final gapcost option for <call> to <xdrop_gap_final>. See blast
   documentation for explanation of option! */
void                gt_blast_process_call_set_xdrop_gap_final(
                                                       GtBlastProcessCall *call,
                                                       double xdrop_gap_final);
/* Sets output format option for <call> to <outfmt>. See blast documentation for
   explanation of option! */
void                gt_blast_process_call_set_outfmt(GtBlastProcessCall *call,
                                                     int outfmt);
/* Sets output format option for <call> to tabular output. See blast
   documentation for explanation of option! */
void                gt_blast_process_call_set_outfmt_tabular(
                                                      GtBlastProcessCall *call);
/* Add string <opt> to options of <call>. */
void                gt_blast_process_call_set_opt(GtBlastProcessCall *call,
                                                  const char *opt);
/* Returns a string representation of the constructed call to blast */
const char*         gt_blast_process_call_get_call(GtBlastProcessCall *call);

/* Run the constructed blast-call and return a pointer to the pipe created with
   popen(3). Caller is responsible for freeing the FILE-pointer with pclose(3).
   Asserts query and db options were set. */
FILE*               gt_blast_process_call_run(GtBlastProcessCall *call,
                                              GtError *err);
/* Free memory of <call>. */
void                gt_blast_process_call_delete(GtBlastProcessCall *call);

#endif
