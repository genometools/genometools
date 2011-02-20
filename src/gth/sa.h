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

#ifndef SA_H
#define SA_H

#include "core/array.h"
#include "core/error.h"
#include "core/file.h"
#include "core/str.h"
#include "core/strand.h"
#include "gth/backtrace_path.h"
#include "gth/bssm_param.h"
#include "gth/gthalphatype.h"
#include "gth/input.h"
#include "gth/stat.h"

/*
  This splice site score is used for protein spliced alignments, because
  for them no splice site scores are computed.
*/

#define UNDEFINED_SPLICE_SITE_SCORE     -1.0

/*
  This struct represents an exon. Thereby, the DNA positions always refer
  to the actual strand.
  Terminology: A DNA position can refer to the forward strand and to the
  actual strand. In the case that the considered position lies on the forward
  strand these are the same. In the case that the considered position lies on
  the reverse complement strand the following relation holds:
  pos_referring_to_actual_strand == total_genomic_length - 1 -
                                    pos_referring_to_forward_strand
*/

typedef struct {
  unsigned long leftgenomicexonborder,    /* the borders of the exons in the */
                rightgenomicexonborder,   /*  genomic sequence               */
                leftreferenceexonborder,  /* the borders of the exons in the */
                rightreferenceexonborder; /* reference sequence              */
  GthDbl exonscore;                       /* = exnscr                        */
} Exoninfo;

typedef struct {
  GthFlt donorsiteprobability,     /* (GS = itrscr) */
         acceptorsiteprobability;  /* (GS = itrscr) */
  GthDbl donorsitescore,           /* (GS = no equivalent) */
         acceptorsitescore;        /* (GS = no equivalent) */
} Introninfo;

/* The following structure bundles all information necessary to represent a
   spliced alignment. Keep in sync with gth_sas_are_equal()! */
typedef struct GthSA GthSA;

/* create a new spliced alignment (sa) */
GthSA*         gth_sa_new(void);
/* create a new sa and set the sequence ids and strand directions */
GthSA*         gth_sa_new_and_set(bool gen_strand_forward,
                                  bool ref_strand_forward,
                                  GthInput *input,
                                  unsigned long gen_file_num,
                                  unsigned long gen_seq_num,
                                  unsigned long ref_file_num,
                                  unsigned long ref_seq_num,
                                  unsigned long call_number,
                                  unsigned long gen_total_length,
                                  unsigned long gen_offset,
                                  unsigned long ref_total_length);
void            gth_sa_set(GthSA*, GthAlphatype ref_alphatype,
                           unsigned long gen_dp_start,
                           unsigned long gen_dp_length);
void            gth_sa_set_gen_dp_length(GthSA*, unsigned long gen_dp_length);
void            gth_sa_delete(GthSA*);
void            gth_sa_show_exons(const GthSA*, GtFile*);
void            gth_sa_get_exons(const GthSA*, GtArray *ranges);
bool            gth_sa_exons_are_forward_and_consecutive(const GthSA*);

/* returns the genomic range of <sa> referring to the forward strand. */
GtRange         gth_sa_range_forward(const GthSA *sa);

/* returns the genomic range of <sa> referring to the actual strand */
GtRange         gth_sa_range_actual(const GthSA *sa);

/* Returns the average splice site probability of spliced alignment <sa>.
   That is, the average probability of all donor and acceptor sites. If the
   spliced alignment contains no introns the average splice site probability is
   zero. */
GthFlt gth_sa_average_splice_site_prob(const GthSA *sa);

/* A spliced alignment <sa> is ``poor'' if one of the following statements
   holds:
   - The number of introns equals zero and the polyA tail of the cDNA/EST is
     at the start of the spliced alignment
   - The number of introns equals one and the average splice site probability
     is less than 1.5 times the minimum average splice site probability
     <minaveragessp>.
   - The number of introns is greater than one and the average splice site
     probability is less than the minimum average splice site probability
     <minaveragessp>, which can be changed by a command line option (see
     manual). */
bool            gth_sa_is_poor(const GthSA *sa, GthFlt minaveragessp);

/* Consider two spliced alignments <saA> and <saB>.
   Spliced alignment <saB> is ``better'' than spliced alignment <saA> if at
   least one of the following statements holds:
   - <saA> has no introns.
   - <saB> has at least one intron and the average splice site probability of
     <saB> is greater than the average splice site probability of <saA>. */
bool            gth_sa_B_is_better_than_A(const GthSA *saA, const GthSA *saB);

/* Returns the left genomic exon border for <exon> of <spliced_alignment>
   (referring to the foward strand). */
unsigned long   gth_sa_left_genomic_exon_border(const GthSA*,
                                                unsigned long exon);

/* Returns the right genomic exon border for <exon> of <spliced_alignment>
   (referring to the foward strand). */
unsigned long   gth_sa_right_genomic_exon_border(const GthSA*,
                                                 unsigned long exon);

double          gth_sa_exon_score(const GthSA*, unsigned long exon);
GtRange         gth_sa_donor_site_range(const GthSA*, unsigned long intron);
GtRange         gth_sa_acceptor_site_range(const GthSA*, unsigned long intron);
/* Return donor site probablity with three positions after decimal point. */
float           gth_sa_donor_site_prob(const GthSA*, unsigned long intron);
/* Return acceptor site probablity with three positions after decimal point. */
float           gth_sa_acceptor_site_prob(const GthSA*,
                                                      unsigned long intron);
unsigned long   gth_sa_genomic_exon_length(const GthSA*, unsigned long exon);
unsigned long   gth_sa_left_intron_border(const GthSA*, unsigned long intron);
unsigned long   gth_sa_right_intron_border(const GthSA*, unsigned long intron);
unsigned long   gth_sa_intron_length(const GthSA *, unsigned long intron);

/* XXX */
GthBacktracePath* gth_sa_backtrace_path(const GthSA*);
Editoperation*  gth_sa_get_editoperations(const GthSA*);
unsigned long   gth_sa_get_editoperations_length(const GthSA*);

unsigned long   gth_sa_indelcount(const GthSA*);
unsigned long   gth_sa_gen_dp_length(const GthSA*);
unsigned long   gth_sa_gen_total_length(const GthSA*);
void            gth_sa_set_gen_total_length(GthSA*, unsigned long);
unsigned long   gth_sa_gen_offset(const GthSA*);
void            gth_sa_set_gen_offset(GthSA*, unsigned long);
unsigned long   gth_sa_ref_total_length(const GthSA*);
void            gth_sa_set_ref_total_length(GthSA*, unsigned long);
unsigned long   gth_sa_gen_dp_start(const GthSA*);
unsigned long   gth_sa_gen_dp_start_show(const GthSA*);
void            gth_sa_set_gen_dp_start(GthSA*, unsigned long);
unsigned long   gth_sa_gen_dp_end(const GthSA*);
unsigned long   gth_sa_gen_dp_end_show(const GthSA*);
unsigned long   gth_sa_gen_file_num(const GthSA*);
void            gth_sa_set_gen_file_num(GthSA*, unsigned long);
unsigned long   gth_sa_gen_seq_num(const GthSA*);
void            gth_sa_set_gen_seq_num(GthSA*, unsigned long);
unsigned long   gth_sa_ref_file_num(const GthSA*);
void            gth_sa_set_ref_file_num(GthSA*, unsigned long);
unsigned long   gth_sa_ref_seq_num(const GthSA*);
void            gth_sa_set_ref_seq_num(GthSA*, unsigned long);
const char*     gth_sa_gen_id(const GthSA*);
GtStr*          gth_sa_gen_id_str(const GthSA*);
void            gth_sa_set_gen_id(GthSA*, const char*);
const char*     gth_sa_ref_id(const GthSA*);
GtStr*          gth_sa_ref_id_str(const GthSA*);
void            gth_sa_set_ref_id(GthSA*, const char*);
GtStr*          gth_sa_gen_md5(const GthSA*);
GtStr*          gth_sa_ref_md5(const GthSA*);
GtStr*          gth_sa_gen_desc(const GthSA*);
GtStr*          gth_sa_ref_desc(const GthSA*);
GtStrand        gth_sa_gen_strand(const GthSA*);
bool            gth_sa_gen_strand_forward(const GthSA*);
char            gth_sa_gen_strand_char(const GthSA*);
void            gth_sa_set_gen_strand(GthSA*, bool forward);
bool            gth_sa_ref_strand_forward(const GthSA*);
char            gth_sa_ref_strand_char(const GthSA*);
void            gth_sa_set_ref_strand(GthSA*, bool forward);
unsigned long   gth_sa_genomiccutoff_start(const GthSA*);
unsigned long   gth_sa_referencecutoff_start(const GthSA*);
unsigned long   gth_sa_eopcutoff_start(const GthSA*);
unsigned long   gth_sa_genomiccutoff_end(const GthSA*);
unsigned long   gth_sa_referencecutoff_end(const GthSA*);
unsigned long   gth_sa_eopcutoff_end(const GthSA*);
void            gth_sa_set_cutoffs_start(GthSA*, Cutoffs*);
void            gth_sa_set_cutoffs_end(GthSA*, Cutoffs*);
GthAlphatype    gth_sa_alphatype(const GthSA*);
const char*     gth_sa_alphastring(const GthSA*);
void            gth_sa_set_alphatype(GthSA*, GthAlphatype);
Exoninfo*       gth_sa_get_exon(const GthSA*, unsigned long);
void            gth_sa_add_exon(GthSA*, Exoninfo*);
unsigned long   gth_sa_num_of_exons(const GthSA*);
Introninfo*     gth_sa_get_intron(const GthSA*, unsigned long);
void            gth_sa_add_intron(GthSA*, Introninfo*);
unsigned long   gth_sa_num_of_introns(const GthSA*);
void            gth_sa_calc_polyAtailpos(GthSA*,
                                         const unsigned char *ref_seq_tran,
                                         GtAlphabet *ref_alpha);
unsigned long   gth_sa_polyAtail_start(const GthSA*);
unsigned long   gth_sa_polyAtail_stop(const GthSA*);
void            gth_sa_set_polyAtail_start(GthSA*, unsigned long);
void            gth_sa_set_polyAtail_stop(GthSA*, unsigned long);
GthFlt          gth_sa_score(const GthSA*);
void            gth_sa_set_score(GthSA*, GthFlt);
GthFlt          gth_sa_coverage(const GthSA*);
void            gth_sa_set_coverage(GthSA*, GthFlt);
bool            gth_sa_genomic_cov_is_highest(const GthSA*);
char            gth_sa_coverage_char(const GthSA*);
void            gth_sa_set_highest_cov(GthSA*, bool genomic);
unsigned long   gth_sa_cumlen_scored_exons(const GthSA*);
void            gth_sa_set_cumlen_scored_exons(GthSA*, unsigned long);
unsigned long   gth_sa_call_number(const GthSA*);
const char*     gth_sa_gff3_target_attribute(GthSA*, bool md5ids);
void            gth_sa_determine_cutoffs(GthSA*, GthCutoffmode leadcutoffsmode,
                                         GthCutoffmode termcutoffsmode,
                                         unsigned long cutoffsminexonlen);
void            gth_sa_cutoff_start(GthSA*);
void            gth_sa_cutoff_end(GthSA*);
void            gth_sa_cutoff_walked_path(GthSA*, const GthPathWalker*,
                                          bool showeops, GtFile*);
void            gth_sa_prepend(GthSA*, const GthBacktracePath*);
void            gth_sa_append(GthSA*, const GthBacktracePath*);
void            gth_sa_remove_zero_base_exons(GthSA*, GthStat *stat);
bool            gth_sa_contains_no_zero_base_exons(const GthSA*);
void            gth_sa_echo_genomic_description(const GthSA*, GthInput*,
                                                GtFile*);
void            gth_sa_echo_reference_description(const GthSA*, GthInput*,
                                                  GtFile*);
void            gth_sa_echo_reference_sequence(const GthSA*, GthInput*,
                                               bool format, GtFile*);
void            gth_sa_echo_alignment(const GthSA *sa,
                                      unsigned long showintronmaxlen,
                                      unsigned long translationtable,
                                      bool wildcardimplosion, GthInput *input,
                                      GtFile *outfp);
/* Fills the alignment lines (<first_line> and <second_line> for DNA alignments
   and <first_line>, <second_line>, and <third_line> for protein alignments). */
unsigned long   gth_sa_get_alignment_lines(const GthSA *sa,
                                           unsigned char **first_line,
                                           unsigned char **second_line,
                                           unsigned char **third_line,
                                           unsigned long translationtable,
                                           GthInput *input);
bool            gth_sa_is_valid(const GthSA*);
void            gth_sa_show(GthSA*, GthInput*, GtFile*);
void            gth_sa_save_ref_md5(GthSA*, GthInput*);

bool            gth_sas_are_equal(const GthSA*, const GthSA*);

#endif
