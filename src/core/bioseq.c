/*
  Copyright (c) 2006-2012 Gordon Gremme <gordon@gremme.org>
  Copyright (c)      2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2006-2012 Center for Bioinformatics, University of Hamburg

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

#include <signal.h>
#include <string.h>
#include "core/bioseq.h"
#include "core/cstr_api.h"
#include "core/disc_distri_api.h"
#include "core/dynalloc.h"
#include "core/encseq.h"
#include "core/fa.h"
#include "core/fasta.h"
#include "core/fasta_reader_fsm.h"
#include "core/fasta_reader_rec.h"
#include "core/fasta_reader_seqit.h"
#include "core/fileutils_api.h"
#include "core/gc_content.h"
#include "core/hashmap_api.h"
#include "core/hashmap-generic.h"
#include "core/ma.h"
#include "core/md5_tab.h"
#include "core/parseutils.h"
#include "core/sig.h"
#include "core/str_array.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "core/xposix.h"

struct GtBioseq {
  bool use_stdin;
  GtStr *sequence_file;
  GtSeq **seqs;
  char **descriptions;
  GtEncseq *encseq;
  GtMD5Tab *md5_tab;
};

/* this global variables are necessary for the signal handler below */
static const char *gt_bioseq_index_filename;

static void remove_indexfile(const char *suffix)
{
  GtStr *fn = gt_str_new_cstr(gt_bioseq_index_filename);
  gt_str_append_cstr(fn, suffix);
  if (gt_file_exists(gt_str_get(fn)))
    gt_xunlink(gt_str_get(fn));
  gt_str_delete(fn);
}

/* removes the incomplete bioseq files */
static void remove_bioseq_files(int sigraised)
{
  remove_indexfile(GT_ENCSEQFILESUFFIX);
  remove_indexfile(GT_DESTABFILESUFFIX);
  remove_indexfile(GT_SSPTABFILESUFFIX);
  remove_indexfile(GT_SDSTABFILESUFFIX);
  remove_indexfile(GT_MD5TABFILESUFFIX);
  remove_indexfile(GT_OISTABFILESUFFIX);
  (void) gt_xsignal(sigraised, SIG_DFL);
  gt_xraise(sigraised);
}

static int construct_bioseq_files(GtBioseq *bs, GtStr *bioseq_indexname,
                                  GtError *err)
{
  GtStr *sequence_filename;
  GtEncseqEncoder *ee;
  GtStrArray *indexfn;
  int had_err = 0;

  gt_error_check(err);

  /* register the signal handler to remove incomplete files upon termination */
  if (!bs->use_stdin) {
    gt_bioseq_index_filename = gt_str_get(bs->sequence_file);
    gt_sig_register_all(remove_bioseq_files);
  }

  /* if stdin is used as input, we need to create a tempfile containing the
     sequence as GtEncseq cannot be built from stdin directly */
  if (bs->use_stdin) {
    GtStr *tmpfilename;
    FILE *tmpfile = NULL;
    int i;
    char buf[BUFSIZ];
    tmpfilename = gt_str_new();
    tmpfile = gt_xtmpfp(tmpfilename);
    gt_assert(tmpfile);
    i = 1;
    while (i > 0) {
      i = fread(buf, 1, BUFSIZ, stdin);
      if (i > 0) fwrite(buf, 1, i, tmpfile);
    }
    gt_fa_xfclose(tmpfile);
    sequence_filename = tmpfilename;
  } else {
    sequence_filename = gt_str_ref(bs->sequence_file);
  }
  gt_assert(gt_str_length(sequence_filename) > 0);
  ee = gt_encseq_encoder_new();
  gt_encseq_encoder_enable_description_support(ee);
  gt_encseq_encoder_enable_md5_support(ee);
  gt_encseq_encoder_enable_multiseq_support(ee);
  gt_encseq_encoder_enable_lossless_support(ee);
  indexfn = gt_str_array_new();
  gt_str_array_add(indexfn, sequence_filename);
  gt_str_delete(sequence_filename);
  had_err = gt_encseq_encoder_encode(ee, indexfn,
                                     gt_str_get(bioseq_indexname), err);
  /* unregister the signal handler */
   if (!bs->use_stdin)
    gt_sig_unregister_all();

  gt_str_array_delete(indexfn);
  gt_encseq_encoder_delete(ee);
  return had_err;
}

static int bioseq_fill(GtBioseq *bs, bool recreate, GtError *err)
{
  GtStr *bioseq_index_file = NULL,
        *bioseq_ois_file = NULL,
        *bioseq_sds_file = NULL,
        *bioseq_md5_file = NULL,
        *bioseq_des_file = NULL;
  int had_err = 0;
  GtStr *bioseq_basename;

  gt_assert(!bs->encseq);

  if (bs->use_stdin)
    bioseq_basename = gt_str_new_cstr("stdin");
  else
    bioseq_basename = bs->sequence_file;

  /* construct file names */
  bioseq_index_file = gt_str_clone(bioseq_basename);
  gt_str_append_cstr(bioseq_index_file, GT_ENCSEQFILESUFFIX);
  bioseq_ois_file = gt_str_clone(bioseq_basename);
  gt_str_append_cstr(bioseq_ois_file, GT_OISTABFILESUFFIX);
  bioseq_sds_file = gt_str_clone(bioseq_basename);
  gt_str_append_cstr(bioseq_sds_file, GT_SDSTABFILESUFFIX);
  bioseq_md5_file = gt_str_clone(bioseq_basename);
  gt_str_append_cstr(bioseq_md5_file, GT_MD5TABFILESUFFIX);
  bioseq_des_file = gt_str_clone(bioseq_basename);
  gt_str_append_cstr(bioseq_des_file, GT_DESTABFILESUFFIX);

  /* construct the bioseq files if necessary */
  if (recreate || bs->use_stdin ||
      !gt_file_exists(gt_str_get(bioseq_index_file)) ||
      !gt_file_exists(gt_str_get(bioseq_ois_file)) ||
      !gt_file_exists(gt_str_get(bioseq_sds_file)) ||
      !gt_file_exists(gt_str_get(bioseq_md5_file)) ||
      !gt_file_exists(gt_str_get(bioseq_des_file)) ||
      gt_file_is_newer(gt_str_get(bs->sequence_file),
                       gt_str_get(bioseq_index_file))) {
    had_err = construct_bioseq_files(bs, bioseq_basename, err);
  }

  if (!had_err) {
    GtEncseqLoader *el = gt_encseq_loader_new();
    gt_encseq_loader_disable_autosupport(el);
    gt_encseq_loader_require_lossless_support(el);
    gt_encseq_loader_require_description_support(el);
    gt_encseq_loader_require_md5_support(el);
    gt_encseq_loader_require_multiseq_support(el);
    bs->encseq = gt_encseq_loader_load(el, gt_str_get(bioseq_basename), err);
    if (bs->encseq == NULL) {
      had_err = -1;
      gt_assert(gt_error_is_set(err));
    }
    gt_encseq_loader_delete(el);
  }
  if (!had_err) {
    gt_assert(bs->encseq);
  }

  /* free */
  if (bs->use_stdin)
    gt_str_delete(bioseq_basename);
  gt_str_delete(bioseq_index_file);
  gt_str_delete(bioseq_ois_file);
  gt_str_delete(bioseq_md5_file);
  gt_str_delete(bioseq_sds_file);
  gt_str_delete(bioseq_des_file);

  return had_err;
}

static GtBioseq* bioseq_new_with_recreate_and_type(GtStr *sequence_file,
                                                   bool recreate, GtError *err)
{
  GtBioseq *bs;
  int had_err = 0;
  gt_error_check(err);
  bs = gt_calloc(1, sizeof *bs);
  if (!strcmp(gt_str_get(sequence_file), "-"))
    bs->use_stdin = true;
  if (!bs->use_stdin && !gt_file_exists(gt_str_get(sequence_file))) {
    gt_error_set(err, "sequence file \"%s\" does not exist or is not readable",
                 gt_str_get(sequence_file));
    had_err = -1;
  }
  if (!had_err) {
    bs->sequence_file = gt_str_ref(sequence_file);
    had_err = bioseq_fill(bs, recreate, err);
  }
  if (had_err) {
    gt_bioseq_delete(bs);
    return NULL;
  }
  gt_assert(bs->encseq);
  bs->descriptions = gt_calloc(gt_encseq_num_of_sequences(bs->encseq),
                               sizeof (char*));
  return bs;
}

GtBioseq* gt_bioseq_new(const char *sequence_file, GtError *err)
{
  GtBioseq *bs;
  GtStr *seqfile;
  gt_error_check(err);
  seqfile = gt_str_new_cstr(sequence_file);
  bs = bioseq_new_with_recreate_and_type(seqfile, false, err);
  gt_str_delete(seqfile);
  return bs;
}

GtBioseq* gt_bioseq_new_recreate(const char *sequence_file, GtError *err)
{
  GtBioseq *bs;
  GtStr *seqfile;
  gt_error_check(err);
  seqfile = gt_str_new_cstr(sequence_file);
  bs = bioseq_new_with_recreate_and_type(seqfile, true, err);
  gt_str_delete(seqfile);
  return bs;
}

GtBioseq* gt_bioseq_new_str(GtStr *sequence_file, GtError *err)
{
  return bioseq_new_with_recreate_and_type(sequence_file, false, err);
}

void gt_bioseq_delete(GtBioseq *bs)
{
  GtUword i;
  if (!bs) return;
  gt_str_delete(bs->sequence_file);
  gt_md5_tab_delete(bs->md5_tab);
  if (bs->descriptions) {
    for (i = 0; i < gt_encseq_num_of_sequences(bs->encseq); i++) {
      gt_free(bs->descriptions[i]);
    }
    gt_free(bs->descriptions);
  }
  gt_encseq_delete(bs->encseq);
  gt_free(bs);
}

GtAlphabet* gt_bioseq_get_alphabet(GtBioseq *bs)
{
  gt_assert(bs);
  return gt_encseq_alphabet(bs->encseq);
}

GtSeq* gt_bioseq_get_seq(GtBioseq *bs, GtUword idx)
{
  GtSeq *seq;
  gt_assert(bs);
  gt_assert(idx < gt_encseq_num_of_sequences(bs->encseq));
  seq = gt_seq_new_own(gt_bioseq_get_sequence(bs, idx),
                       gt_bioseq_get_sequence_length(bs, idx),
                       gt_encseq_alphabet(bs->encseq));
  gt_seq_set_description(seq, gt_bioseq_get_description(bs, idx));
  return seq;
}

GtSeq* gt_bioseq_get_seq_range(GtBioseq *bs, GtUword idx,
                               GtUword start, GtUword end)
{
  GtSeq *seq;
  gt_assert(bs);
  gt_assert(idx < gt_encseq_num_of_sequences(bs->encseq));
  gt_assert(end >= start);
  gt_assert(end - start + 1 > gt_encseq_seqlength(bs->encseq, idx));
  seq = gt_seq_new_own(gt_bioseq_get_sequence_range(bs, idx, start, end),
                       end - start + 1,
                       gt_encseq_alphabet(bs->encseq));
  gt_seq_set_description(seq, gt_bioseq_get_description(bs, idx));
  return seq;
}

const char* gt_bioseq_get_description(GtBioseq *bs, GtUword idx)
{
  const char *desc;
  char *mydesc;
  GtUword desclen;
  gt_assert(bs && bs->encseq);
  gt_assert(idx < gt_encseq_num_of_sequences(bs->encseq));
  if (!(mydesc = bs->descriptions[idx])) {
    desc = gt_encseq_description(bs->encseq, &desclen, idx);
    mydesc = gt_calloc(desclen + 1, sizeof (char));
    strncpy(mydesc, desc, desclen);
    bs->descriptions[idx] = mydesc;
  }
  return (const char*) mydesc;
}

char gt_bioseq_get_char(const GtBioseq *bs, GtUword index,
                        GtUword position)
{
  GtUword startpos;
  gt_assert(bs);
  gt_assert(index < gt_encseq_num_of_sequences(bs->encseq));
  startpos = gt_encseq_seqstartpos(bs->encseq, index);
  return gt_encseq_get_decoded_char(bs->encseq, startpos + position,
                                    GT_READMODE_FORWARD);
}

char* gt_bioseq_get_sequence(const GtBioseq *bs, GtUword idx)
{
  char *out;
  GtUword startpos;
  gt_assert(bs);
  gt_assert(idx < gt_encseq_num_of_sequences(bs->encseq));
  out = gt_calloc(gt_encseq_seqlength(bs->encseq, idx), sizeof (char));
  startpos = gt_encseq_seqstartpos(bs->encseq, idx);
  gt_encseq_extract_decoded(bs->encseq, out, startpos,
                            startpos
                              + gt_encseq_seqlength(bs->encseq, idx) - 1);
  return out;
}

char* gt_bioseq_get_sequence_range(const GtBioseq *bs, GtUword idx,
                                   GtUword start, GtUword end)
{
  char *out;
  GtUword startpos;
  gt_assert(bs);
  gt_assert(idx < gt_encseq_num_of_sequences(bs->encseq) && end >= start);
  out = gt_malloc((end - start + 1) * sizeof (char));
  startpos = gt_encseq_seqstartpos(bs->encseq, idx);
  gt_encseq_extract_decoded(bs->encseq, out, startpos + start, startpos + end);
  return out;
}

GtUchar gt_bioseq_get_encoded_char(const GtBioseq *bs, GtUword index,
                                   GtUword position)
{
  GtUword startpos;
  gt_assert(bs);
  gt_assert(index < gt_encseq_num_of_sequences(bs->encseq));
  startpos = gt_encseq_seqstartpos(bs->encseq, index);
  return gt_encseq_get_encoded_char(bs->encseq, startpos + position,
                                    GT_READMODE_FORWARD);
}

void gt_bioseq_get_encoded_sequence(const GtBioseq *bs, GtUchar *out,
                                    GtUword idx)
{
  GtUword startpos;
  gt_assert(bs);
  gt_assert(idx < gt_encseq_num_of_sequences(bs->encseq));
  startpos = gt_encseq_seqstartpos(bs->encseq, idx);
  gt_encseq_extract_encoded(bs->encseq, out, startpos,
                            startpos
                              + gt_encseq_seqlength(bs->encseq, idx) - 1);
}

void gt_bioseq_get_encoded_sequence_range(const GtBioseq *bs, GtUchar *out,
                                          GtUword idx,
                                          GtUword start,
                                          GtUword end)
{
  GtUword startpos;
  gt_assert(bs);
  gt_assert(idx < gt_encseq_num_of_sequences(bs->encseq) && end >= start);
  startpos = gt_encseq_seqstartpos(bs->encseq, idx);
  gt_encseq_extract_encoded(bs->encseq, out, startpos + start, startpos + end);
}

const char* gt_bioseq_get_md5_fingerprint(GtBioseq *bs, GtUword idx)
{
  gt_assert(bs && idx < gt_bioseq_number_of_sequences(bs));
  if (!bs->md5_tab) {
    bs->md5_tab = gt_encseq_get_md5_tab(bs->encseq, NULL);
  }
  gt_assert(gt_md5_tab_get(bs->md5_tab, idx));
  return gt_md5_tab_get(bs->md5_tab, idx);
}

const char* gt_bioseq_filename(const GtBioseq *bs)
{
  gt_assert(bs);
  return gt_str_get(bs->sequence_file);
}

GtUword gt_bioseq_get_sequence_length(const GtBioseq *bs,
                                            GtUword idx)
{
  gt_assert(bs);
  return gt_encseq_seqlength(bs->encseq, idx);
}

GtUword gt_bioseq_get_total_length(const GtBioseq *bs)
{
  gt_assert(bs);
  return gt_encseq_total_length(bs->encseq)
           - gt_encseq_num_of_sequences(bs->encseq) + 1;
}

GtUword gt_bioseq_number_of_sequences(GtBioseq *bs)
{
  gt_assert(bs);
  return gt_encseq_num_of_sequences(bs->encseq);
}

GtUword gt_bioseq_md5_to_index(GtBioseq *bs, const char *md5)
{
  gt_assert(bs && md5 && gt_encseq_has_md5_support(bs->encseq));
  if (!bs->md5_tab) {
    bs->md5_tab = gt_encseq_get_md5_tab(bs->encseq, NULL);
  }
  return gt_md5_tab_map(bs->md5_tab, md5);
}

void gt_bioseq_show_as_fasta(GtBioseq *bs, GtUword width, GtFile *outfp)
{
  GtUword i;

  gt_assert(bs);

  for (i = 0; i < gt_bioseq_number_of_sequences(bs); i++) {
    char *seq = gt_bioseq_get_sequence(bs, i);
    gt_fasta_show_entry(gt_bioseq_get_description(bs, i),
                        seq,
                        gt_bioseq_get_sequence_length(bs, i), width, outfp);
    gt_free(seq);
  }
}

void gt_bioseq_show_sequence_as_fasta(GtBioseq *bs, GtUword seqnum,
                                      GtUword width, GtFile *outfp)
{
  char *seq = NULL;
  gt_assert(bs);
  gt_assert(seqnum < gt_bioseq_number_of_sequences(bs));
  seq = gt_bioseq_get_sequence(bs, seqnum);

  gt_fasta_show_entry(gt_bioseq_get_description(bs, seqnum),
                      seq,
                      gt_bioseq_get_sequence_length(bs, seqnum), width, outfp);

  gt_free(seq);
}

void gt_bioseq_show_gc_content(GtBioseq *bs, GtFile *outfp)
{
  gt_assert(bs);
  if (gt_alphabet_is_dna(gt_encseq_alphabet(bs->encseq))) {
    GtUword i, GT_UNUSED purecharlen;
    GtStr *str = gt_str_new();
    purecharlen = gt_encseq_total_length(bs->encseq)
                    - gt_encseq_num_of_sequences(bs->encseq) + 1;
    for (i=0; i < gt_encseq_num_of_sequences(bs->encseq); i++) {
      char *tmp;
      tmp = gt_bioseq_get_sequence(bs, i);
      gt_str_append_cstr(str, tmp);
      gt_free(tmp);
    }
    gt_assert(gt_str_length(str) == purecharlen);
    gt_file_xprintf(outfp, "showing GC-content for sequence file \"%s\"\n",
                    gt_str_get(bs->sequence_file));
    gt_gc_content_show(gt_str_get(str),
                       gt_str_length(str),
                       gt_encseq_alphabet(bs->encseq),
                       outfp);
    gt_str_delete(str);
  }
}

void gt_bioseq_show_stat(GtBioseq *bs, GtFile *outfp)
{
  GtUword i, num_of_seqs;
  gt_assert(bs);
  num_of_seqs = gt_bioseq_number_of_sequences(bs);
  gt_file_xprintf(outfp, "showing statistics for sequence file \"%s\"\n",
                  gt_str_get(bs->sequence_file));
  gt_file_xprintf(outfp, "number of sequences: "GT_WU"\n", num_of_seqs);
  gt_file_xprintf(outfp, "total length: "GT_WU"\n",
                    gt_encseq_total_length(bs->encseq)
                      - gt_encseq_num_of_sequences(bs->encseq) + 1);
  for (i = 0; i < num_of_seqs; i++) {
    gt_file_xprintf(outfp, "sequence #"GT_WU" length: "GT_WU"\n", i+1,
                    gt_bioseq_get_sequence_length(bs, i));
  }
}

void gt_bioseq_show_seqlengthdistri(GtBioseq *bs, GtFile *outfp)
{
  GtDiscDistri *d;
  GtUword i;
  gt_assert(bs);
  d = gt_disc_distri_new();
  for (i = 0; i < gt_bioseq_number_of_sequences(bs); i++)
    gt_disc_distri_add(d, gt_bioseq_get_sequence_length(bs, i));
  gt_file_xprintf(outfp, "sequence length distribution:\n");
  gt_disc_distri_show(d, outfp);
  gt_disc_distri_delete(d);
}
