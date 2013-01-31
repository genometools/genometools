/*
  Copyright (c) 2003-2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef INPUT_H
#define INPUT_H

#include "core/array.h"
#include "core/range.h"
#include "core/score_matrix.h"
#include "core/str_array.h"
#include "gth/duplicate_check.h"
#include "gth/gthoutput.h"
#include "gth/gthalphatype.h"
#include "gth/seq_con.h"

typedef struct GthInput GthInput;

typedef int (*GthInputFilePreprocessor)(GthInput *input,
                                        bool gthconsensus,
                                        bool noautoindex,
                                        bool skipindexcheck,
                                        bool maskpolyAtails,
                                        bool online,
                                        bool inverse,
                                        const char *progname,
                                        unsigned int translationtable,
                                        GthOutput *out, GtError*);

GthInput*      gth_input_new(GthInputFilePreprocessor file_preprocessor,
                             GthSeqConConstructor seq_con_constructor);
int            gth_input_preprocess(GthInput*,
                                    bool gthconsensus,
                                    bool noautoindex,
                                    bool createindicesonly,
                                    bool skipindexcheck,
                                    bool maskpolyAtails,
                                    bool online,
                                    bool inverse,
                                    const char *progname,
                                    char *scorematrixfile,
                                    unsigned int translationtable,
                                    GthDuplicateCheck duplicate_check,
                                    GthOutput *out, GtError*);
void           gth_input_add_genomic_file(GthInput*, const char *filename);
void           gth_input_add_cdna_file(GthInput*, const char *filename);
void           gth_input_add_protein_file(GthInput*, const char *filename);
void           gth_input_add_reference_file(GthInput*, const char *filename,
                                            GthAlphatype);
const char*    gth_input_get_genomic_filename(const GthInput*,
                                              unsigned long gen_file_num);
const char*    gth_input_get_reference_filename(const GthInput*,
                                                unsigned long ref_file_num);
GthAlphatype   gth_input_get_alphatype(const GthInput*, unsigned long);
bool           gth_input_ref_file_is_dna(const GthInput*, unsigned long);
bool           gth_input_md5ids(const GthInput*);
const
unsigned char* gth_input_original_genomic_sequence(GthInput*,
                                                   unsigned long filenum,
                                                   bool forward);
void           gth_input_echo_genomic_description(GthInput*,
                                                  unsigned long filenum,
                                                  unsigned long seqnum,
                                                  GtFile*);
void           gth_input_echo_reference_description(GthInput*,
                                                    unsigned long filenum,
                                                    unsigned long seqnum,
                                                    GtFile*);
void           gth_input_echo_reference_sequence(GthInput*, bool format,
                                                 unsigned long filenum,
                                                 unsigned long seqnum,
                                                 bool forward,
                                                 GtFile*);
void           gth_input_get_genomic_description(GthInput*, GtStr *desc,
                                                 unsigned long filenum,
                                                 unsigned long seqnum);
void           gth_input_save_gen_id(GthInput*, GtStr *id,
                                     unsigned long file_num,
                                     unsigned long seq_num);
/* save the genomic identifier used for output */
void           gth_input_save_gen_identifier(GthInput*, GtStr *id,
                                             unsigned long file_num,
                                             unsigned long seq_num);
void           gth_input_save_ref_id(GthInput*, GtStr *id,
                                     unsigned long file_num,
                                     unsigned long seq_num);
void           gth_input_save_gen_md5(GthInput*, GtStr **md5,
                                      unsigned long file_num,
                                      unsigned long seq_num);
void           gth_input_save_ref_md5(GthInput*, GtStr **md5,
                                      unsigned long file_num,
                                      unsigned long seq_num);
void           gth_input_save_gen_desc(GthInput*, GtStr **desc,
                                       unsigned long file_num,
                                       unsigned long seq_num);
void           gth_input_save_ref_desc(GthInput*, GtStr **desc,
                                       unsigned long file_num,
                                       unsigned long seq_num);
unsigned long  gth_input_num_of_gen_files(const GthInput*);
unsigned long  gth_input_num_of_ref_files(const GthInput*);
unsigned long  gth_input_num_of_gen_seqs(GthInput*, unsigned long filenum);
unsigned long  gth_input_num_of_ref_seqs(GthInput*, unsigned long filenum);
unsigned long  gth_input_genomic_file_total_length(GthInput*,
                                                   unsigned long filenum);
GtRange        gth_input_get_relative_genomic_range(GthInput*,
                                                    unsigned long filenum,
                                                    unsigned long seqnum);
GtRange        gth_input_get_genomic_range(GthInput*,
                                           unsigned long filenum,
                                           unsigned long seqnum);
GtRange        gth_input_get_reference_range(GthInput*,
                                             unsigned long filenum,
                                             unsigned long seqnum);

#define        gth_input_load_genomic_file(input, gen_file_num, translate) \
               gth_input_load_genomic_file_func(input, gen_file_num, translate,\
                                                __FILE__, __LINE__)
void           gth_input_load_genomic_file_func(GthInput *input,
                                                unsigned long gen_file_num,
                                                bool translate,
                                                const char *src_file,
                                                int src_line);
#define        gth_input_load_reference_file(input, ref_file_num, translate) \
               gth_input_load_reference_file_func(input, ref_file_num,\
                                                  translate, __FILE__, __LINE__)
void           gth_input_load_reference_file_func(GthInput *input,
                                                  unsigned long ref_file_num,
                                                  bool translate,
                                                  const char *src_file,
                                                  int src_line);
int            gth_input_load_scorematrix(GthInput*, char *scorematrixfile,
                                          GthOutput *out, GtError *err);
GtStr*         gth_input_proteinsmap(const GthInput*);
GtStr*         gth_input_bssmfile(const GthInput*);
const char*    gth_input_bssmfilename(const GthInput*);
GtScoreMatrix* gth_input_score_matrix(const GthInput*);
GtAlphabet*    gth_input_score_matrix_alpha(const GthInput*);
void           gth_input_set_forward_only(GthInput*);
void           gth_input_set_reverse_only(GthInput*);
bool           gth_input_forward(const GthInput*);
bool           gth_input_reverse(const GthInput*);
bool           gth_input_both(const GthInput*);
GthAlphatype   gth_input_overall_alphatype(const GthInput*);
long           gth_input_determine_genomic_file_index(const GthInput*,
                                                      const char *filename);
long           gth_input_determine_reference_file_index(const GthInput*,
                                                        const char *filename);
int            gth_input_set_and_check_substring_spec(GthInput*, GtError*);
bool           gth_input_use_substring_spec(const GthInput*);
unsigned long  gth_input_genomic_substring_from(const GthInput*);
unsigned long  gth_input_genomic_substring_to(const GthInput*);
/* for option parsing */
unsigned long* gth_input_genomicfrompos_ptr(GthInput*);
/* for option parsing */
unsigned long* gth_input_genomicwidth_ptr(GthInput*);
/* for option parsing */
unsigned long* gth_input_genomictopos_ptr(GthInput*);

int            gth_input_make_indices(GthInput*, const char *progname,
                                      GtError*);

/* The following function frees the current genomic and reference virtual tree
   of the Gthinptinfo structure, but not the BSSM paramter file. */
void           gth_input_delete_current(GthInput*);

/* The following function frees the current genomic and reference virtual tree
   of the GthInput structure and the BSSM paramter file. */
void           gth_input_delete_complete(GthInput*);

/* deprecated interfaces */
GthSeqCon*           gth_input_current_gen_seq_con(GthInput*);
GthSeqCon*           gth_input_current_ref_seq_con(GthInput*);
const unsigned char* gth_input_current_gen_seq_tran(const GthInput*);
const unsigned char* gth_input_current_gen_seq_tran_rc(const GthInput*);
const unsigned char* gth_input_current_gen_seq_orig(const GthInput*);
const unsigned char* gth_input_current_gen_seq_orig_rc(const GthInput*);
const unsigned char* gth_input_current_ref_seq_tran(const GthInput*);
const unsigned char* gth_input_current_ref_seq_tran_rc(const GthInput*);
const unsigned char* gth_input_current_ref_seq_orig(const GthInput*);
const unsigned char* gth_input_current_ref_seq_orig_rc(const GthInput*);
GtAlphabet*          gth_input_current_gen_alphabet(GthInput*);
GtAlphabet*          gth_input_current_ref_alphabet(GthInput*);

#endif
