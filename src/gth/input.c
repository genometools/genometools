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

#include <ctype.h>
#include <math.h>
#include "core/fileutils_api.h"
#include "core/md5_seqid.h"
#include "core/safearith.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/regular_seqid.h"
#include "gth/default.h"
#include "gth/desc_cache.h"
#include "gth/gthdef.h"
#include "gth/input.h"
#include "gth/md5_cache.h"
#include "gth/parse_options.h"

#define GTHFORWARD  1        /* run program to match forward */
#define GTHREVERSE  (1 << 1) /* run program to match reverse */

#define UPPERBLASTMATRIXALPHABET "ARNDCQEGHILKMFPSTWYVBZX*"

struct GthInput {
  GthInputFilePreprocessor file_preprocessor;
  GtStrArray *genomicfiles,       /* pointer to genomic file names */
             *referencefiles;     /* pointer to reference file names */
  GtArray *alphatypes;            /* stores the alphabet types of the reference
                                     files */
  GthAlphatype overall_alphatype; /* the overall alphabet type (dna, protein, or
                                     mixed) */
  unsigned long gen_file_num,     /* number of genomic file which is currently
                                     mapped into memory */
                ref_file_num;     /* number of reference file which is currently
                                     mapped into memory */
  GthSeqConConstructor seq_con_constructor;
  GthSeqCon *gen_seq_con;         /* containing the preprocessed genomic
                                     sequence */
  GthSeqCon *ref_seq_con;         /* containing the reference sequence */
  GtStr *proteinsmap,             /* name of protein symbol mapping */
        *bssmfile;                /* file for bssm parameter */
  GtScoreMatrix *score_matrix;    /* contains the amino acid substitution matrix
                                     (if necessary) */
  GtAlphabet *score_matrix_alpha; /* the alphabet used for the scoring matrix */

  unsigned int searchmode;        /* stores bits GTHFORWARD and GTHREVERSE */
  bool use_substring_spec,
       genomic_translate,
       reference_translate;
  unsigned long genomicfrompos,   /* analyse genomic seq. from this position */
                genomicwidth,     /* analyse only this width of genomic seq. */
                genomictopos;     /* = genomicfrompos + genomicwidth - 1 */
  bool md5ids,
       use_md5_cache,
       use_desc_cache;
  GthMD5Cache *gen_md5_cache,
              *ref_md5_cache;
  GthDescCache *gen_desc_cache,
               *ref_desc_cache;
};

GthInput *gth_input_new(GthInputFilePreprocessor file_preprocessor,
                        GthSeqConConstructor seq_con_constructor)
{
  GthInput *input = gt_calloc(1, sizeof *input);
  gt_assert(seq_con_constructor);
  input->file_preprocessor = file_preprocessor;
  input->genomicfiles = gt_str_array_new();
  input->referencefiles = gt_str_array_new();
  input->alphatypes = gt_array_new(sizeof (GthAlphatype));
  input->overall_alphatype = UNDEF_ALPHA;
  input->gen_file_num = GT_UNDEF_ULONG;
  input->ref_file_num = GT_UNDEF_ULONG;
  input->seq_con_constructor = seq_con_constructor;
  input->proteinsmap = gt_str_new_cstr(GTH_DEFAULT_PROTEINSMAP);
  input->bssmfile = gt_str_new();
  input->searchmode = GTHFORWARD | GTHREVERSE;
  input->genomic_translate = GT_UNDEF_BOOL;
  input->reference_translate = GT_UNDEF_BOOL;
  return input;
}

static void create_md5_cache_files(GthInput *input)
{
  GthMD5Cache *md5_cache;
  const char *filename;
  GthSeqCon *seq_con;
  GtStr *indexname;
  unsigned long i;
  gt_assert(input);
  indexname = gt_str_new();
  for (i = 0; i < gt_str_array_size(input->genomicfiles); i++) {
    filename = gth_input_get_genomic_filename(input, i);
    gt_str_set(indexname, filename);
    gt_str_append_char(indexname, '.');
    gt_str_append_cstr(indexname, DNASUFFIX);
    seq_con = input->seq_con_constructor(gt_str_get(indexname),
                                         false, true,false);
    md5_cache = gth_md5_cache_new(filename, seq_con);
    gth_md5_cache_delete(md5_cache);
    gth_seq_con_delete(seq_con);
  }
  for (i = 0; i < gt_str_array_size(input->referencefiles); i++) {
    GthAlphatype alphatype = gth_input_get_alphatype(input, i);
    filename = gth_input_get_reference_filename(input, i);
    gt_str_set(indexname, filename);
    gt_str_append_char(indexname, '.');
    gt_str_append_cstr(indexname, alphatype == DNA_ALPHA
                                  ? DNASUFFIX
                                  : gt_str_get(input->proteinsmap));
    seq_con = input->seq_con_constructor(gt_str_get(indexname),
                                         false, true, false);
    md5_cache = gth_md5_cache_new(filename, seq_con);
    gth_md5_cache_delete(md5_cache);
    gth_seq_con_delete(seq_con);
  }
  gt_str_delete(indexname);
}

int gth_input_preprocess(GthInput *input,
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
                         GthOutput *out, GtError *err)
{
  int had_err;
  gt_error_check(err);
  gt_assert(input);
  had_err = input->file_preprocessor(input, gthconsensus, noautoindex,
                                     skipindexcheck, maskpolyAtails, online,
                                     inverse, progname, translationtable, out,
                                     err);
  if (!had_err)
    had_err = gth_input_load_scorematrix(input, scorematrixfile, out, err);
  if (!had_err) {
    input->md5ids = out->md5ids;
    input->use_md5_cache = out->md5ids || duplicate_check == GTH_DC_SEQ ||
                           duplicate_check == GTH_DC_BOTH;
    input->use_desc_cache = duplicate_check == GTH_DC_DESC ||
                            duplicate_check == GTH_DC_BOTH;
  }
  if (!had_err && input->use_md5_cache && createindicesonly) {
    /* for performance reasons (mapping of all index files), this is only done
       if <createindicesonly> is <true>. otherwise, this is done automatically,
       if the corresponding cache is accessed. */
    create_md5_cache_files(input);
  }
  return had_err;
}

void gth_input_add_genomic_file(GthInput *input, const char *filename)
{
  gt_assert(input && filename);
  gt_str_array_add_cstr(input->genomicfiles, filename);
}

void gth_input_add_cdna_file(GthInput *input, const char *filename)
{
  gt_assert(input && filename);
  gth_input_add_reference_file(input, filename, DNA_ALPHA);
}

void gth_input_add_protein_file(GthInput *input, const char *filename)
{
  gt_assert(input && filename);
  gth_input_add_reference_file(input, filename, PROTEIN_ALPHA);
}

void gth_input_add_reference_file(GthInput *input, const char *filename,
                                  GthAlphatype alphatype)
{
  gt_assert(input && filename);
  gt_assert(alphatype == DNA_ALPHA || alphatype == PROTEIN_ALPHA);
  gt_str_array_add_cstr(input->referencefiles, filename);
  gt_array_add(input->alphatypes, alphatype);
  /* update overall alphatype */
  if (input->overall_alphatype == UNDEF_ALPHA)
    input->overall_alphatype = alphatype;
  else if (input->overall_alphatype != alphatype)
    input->overall_alphatype = MIXED_ALPHA;
}

const char* gth_input_get_genomic_filename(const GthInput *input,
                                           unsigned long gen_file_num)
{
  gt_assert(input);
  return gt_str_array_get(input->genomicfiles, gen_file_num);
}

const char* gth_input_get_reference_filename(const GthInput *input,
                                             unsigned long ref_file_num)
{
  gt_assert(input);
  return gt_str_array_get(input->referencefiles, ref_file_num);
}

GthAlphatype gth_input_get_alphatype(const GthInput *input,
                                     unsigned long ref_file_num)
{
  gt_assert(input);
  return *(GthAlphatype*) gt_array_get(input->alphatypes, ref_file_num);
}

bool gth_input_ref_file_is_dna(const GthInput *input,
                               unsigned long ref_file_num)
{
  gt_assert(input);
  if (gth_input_get_alphatype(input, ref_file_num) == DNA_ALPHA)
    return true;
  return false;
}

bool gth_input_md5ids(const GthInput *input)
{
  gt_assert(input);
  return input->md5ids;
}

const unsigned char* gth_input_original_genomic_sequence(GthInput *input,
                                                         GT_UNUSED
                                                         unsigned long filenum,
                                                         bool forward)
{
  gt_assert(input);
  gt_assert(input->gen_file_num == filenum);
  if (forward)
    return gth_seq_con_get_orig_seq(input->gen_seq_con, 0);
  else {
    gt_assert(input->searchmode & GTHREVERSE);
    return gth_seq_con_get_orig_seq_rc(input->gen_seq_con, 0);
  }
}

void gth_input_echo_genomic_description(GthInput *input,
                                        GT_UNUSED unsigned long filenum,
                                        unsigned long seqnum, GtFile *outfp)
{
  gt_assert(input);
  gt_assert(input->gen_file_num == filenum);
  gth_seq_con_echo_description(input->gen_seq_con, seqnum, outfp);
}

void gth_input_echo_reference_description(GthInput *input,
                                          GT_UNUSED unsigned long filenum,
                                          unsigned long seqnum,
                                          GtFile *outfp)
{
  gt_assert(input);
  gt_assert(input->ref_file_num == filenum);
  gth_seq_con_echo_description(input->ref_seq_con, seqnum, outfp);
}

static void format_reference_seq(unsigned char *seq, unsigned long len,
                                 GtFile *outfp)
{
  unsigned long i, j, tennercount;
  int width = 1;
  bool showcharnum = true;

  gt_assert(seq && len);

  /* determine necessary width for character numbering */
  if (ALIGNMENTLINEWIDTH < len)
    width = ceil(log10(len - ALIGNMENTLINEWIDTH));

  /* show sequence */
  for (i = 0, j = 0, tennercount = 0; /* nothing */; i++) {
    if (showcharnum) {
      gt_file_xprintf(outfp, "  %*lu  ", width, i + OUTPUTOFFSET);
      showcharnum = false;
    }
    gt_file_xfputc(toupper(seq[i]), outfp);
    if (i == len - 1) {
      gt_file_xfputc('\n', outfp);
      break;
    }
    j++;
    if (j >= ALIGNMENTLINEWIDTH) {
      gt_file_xfputc('\n', outfp);
      j = 0;
      tennercount = 0;
      showcharnum = true;
    }
    else {
      tennercount++;
      if (tennercount == 10) {
        gt_file_xfputc(' ', outfp);
        tennercount = 0;
      }
    }
  }

  gt_file_xfputc('\n', outfp);
}

void gth_input_echo_reference_sequence(GthInput *input, bool format,
                                       GT_UNUSED unsigned long filenum,
                                       unsigned long seqnum, bool forward,
                                       GtFile *outfp)
{
  unsigned char *refseq;
  unsigned long i, reflength;
  gt_assert(input);
  gt_assert(input->ref_file_num == filenum);

  /* get reference sequence */
  if (forward)
    refseq = gth_seq_con_get_orig_seq(input->ref_seq_con, seqnum);
  else
    refseq = gth_seq_con_get_orig_seq_rc(input->ref_seq_con, seqnum);

  /* output reference sequence */
  reflength = gth_seq_con_get_length(input->ref_seq_con, seqnum);
  if (format)
    format_reference_seq(refseq, reflength, outfp);
  else {
    for (i = 0; i < reflength; i++)
      gt_file_xfputc(refseq[i], outfp);
  }
}

void gth_input_get_genomic_description(GthInput *input, GtStr *description,
                                       GT_UNUSED unsigned long filenum,
                                       unsigned long seqnum)
{
  gt_assert(input && description);
  gt_assert(input->gen_file_num == filenum);
  gth_seq_con_get_description(input->gen_seq_con, seqnum, description);
}

static void save_sequenceid(GtStr *sequenceid, GthSeqCon *seqcol,
                            unsigned long seqnum)
{
  GtStr *description;

  /* sequence number is defined */
  gt_assert(seqnum != GT_UNDEF_ULONG);

  description = gt_str_new();
  gth_seq_con_get_description(seqcol, seqnum, description);

  gt_regular_seqid_save(sequenceid, description);

  gt_str_delete(description);
}

void gth_input_save_gen_id(GthInput *input, GtStr *id,
                           GT_UNUSED unsigned long file_num,
                           unsigned long seq_num)
{
  gt_assert(input && id);
  gt_assert(input->gen_file_num == file_num);
  save_sequenceid(id, input->gen_seq_con, seq_num);
}

void gth_input_save_gen_identifier(GthInput *input, GtStr *id,
                                   GT_UNUSED unsigned long file_num,
                                   unsigned long seq_num)
{
  gt_assert(input && id);
  gt_assert(input->gen_file_num == file_num);
  if (input->md5ids) {
    GtStr *md5;
    gt_str_append_cstr(id, GT_MD5_SEQID_PREFIX);
    md5 = gth_md5_cache_get(input->gen_md5_cache, seq_num);
    gt_str_append_str(id, md5);
    gt_str_delete(md5);
    gt_str_append_char(id, ':');
  }
  save_sequenceid(id, input->gen_seq_con, seq_num);
}

void gth_input_save_ref_id(GthInput *input, GtStr *id,
                           GT_UNUSED unsigned long file_num,
                           unsigned long seq_num)
{
  gt_assert(input && id);
  gt_assert(input->ref_file_num == file_num);
  save_sequenceid(id, input->ref_seq_con, seq_num);
}

void gth_input_save_gen_md5(GthInput *input, GtStr **md5,
                            GT_UNUSED unsigned long file_num,
                            unsigned long seq_num)
{
  gt_assert(input && input->gen_file_num == file_num);
  gt_assert(md5 && !*md5);
  if (input->use_md5_cache)
    *md5 = gth_md5_cache_get(input->gen_md5_cache, seq_num);
}

void gth_input_save_ref_md5(GthInput *input, GtStr **md5,
                            GT_UNUSED unsigned long file_num,
                            unsigned long seq_num)
{
  gt_assert(input && input->ref_file_num == file_num);
  gt_assert(md5 && !*md5);
  if (input->use_md5_cache)
    *md5 = gth_md5_cache_get(input->ref_md5_cache, seq_num);
}

void gth_input_save_gen_desc(GthInput *input, GtStr **desc,
                             GT_UNUSED unsigned long file_num,
                             unsigned long seq_num)
{
  gt_assert(input && input->gen_file_num == file_num);
  gt_assert(desc && !*desc);
  if (input->use_desc_cache)
    *desc = gth_desc_cache_get(input->gen_desc_cache, seq_num);
}

void gth_input_save_ref_desc(GthInput *input, GtStr **desc,
                             GT_UNUSED unsigned long file_num,
                             unsigned long seq_num)
{
  gt_assert(input && input->ref_file_num == file_num);
  gt_assert(desc && !*desc);
  if (input->use_desc_cache)
    *desc = gth_desc_cache_get(input->ref_desc_cache, seq_num);
}

unsigned long gth_input_num_of_gen_files(const GthInput *input)
{
  gt_assert(input);
  return gt_str_array_size(input->genomicfiles);
}

unsigned long gth_input_num_of_ref_files(const GthInput *input)
{
  gt_assert(input);
  return gt_str_array_size(input->referencefiles);
}

unsigned long gth_input_genomic_file_total_length(GthInput *input,
                                                  GT_UNUSED
                                                  unsigned long filenum)
{
  gt_assert(input);
  gt_assert(input->gen_file_num == filenum);
  return gth_seq_con_total_length(input->gen_seq_con);
}

unsigned long gth_input_num_of_gen_seqs(GthInput *input,
                                        GT_UNUSED unsigned long filenum)
{
  gt_assert(input);
  gt_assert(input->gen_file_num == filenum);
  return gth_seq_con_num_of_seqs(input->gen_seq_con);
}

unsigned long gth_input_num_of_ref_seqs(GthInput *input,
                                        GT_UNUSED unsigned long filenum)
{
  gt_assert(input);
  gt_assert(input->ref_file_num == filenum);
  return gth_seq_con_num_of_seqs(input->ref_seq_con);
}

GtRange gth_input_get_relative_genomic_range(GthInput *input,
                                             GT_UNUSED unsigned long filenum,
                                             unsigned long seqnum)
{
  gt_assert(input);
  gt_assert(input->gen_file_num == filenum);
  return gth_seq_con_get_relative_range(input->gen_seq_con, seqnum);
}

GtRange gth_input_get_genomic_range(GthInput *input,
                                    GT_UNUSED unsigned long filenum,
                                    unsigned long seqnum)
{
  gt_assert(input);
  gt_assert(input->gen_file_num == filenum);
  return gth_seq_con_get_range(input->gen_seq_con, seqnum);
}

GtRange gth_input_get_reference_range(GthInput *input,
                                      GT_UNUSED unsigned long filenum,
                                      unsigned long seqnum)
{
  gt_assert(input);
  gt_assert(input->ref_file_num == filenum);
  return gth_seq_con_get_range(input->ref_seq_con, seqnum);
}

#define INPUT_DEBUG 0

void gth_input_load_genomic_file_func(GthInput *input,
                                      unsigned long gen_file_num,
                                      bool translate,
                                      GT_UNUSED const char *src_file,
                                      GT_UNUSED int src_line)
{
  /* valid genomic file number */
  gt_assert(input && gen_file_num < gt_str_array_size(input->genomicfiles));

#if INPUT_DEBUG
  printf("load genomic   (file %s, line %d): gen_file_num: %lu, translate=%s\n",
         src_file, src_line, gen_file_num, translate ? "true" : "false");
#endif

  if (input->gen_file_num != gen_file_num) {
    const char *genomic_filename;
    GtStr *indexname;

    /* free old genomic file */
    if (input->gen_file_num != GT_UNDEF_ULONG) {
      /* in this case a sequence collection has been loaded -> free it */
      gth_seq_con_delete(input->gen_seq_con);
      /* delete caches */
      gth_md5_cache_delete(input->gen_md5_cache);
      gth_desc_cache_delete(input->gen_desc_cache);
    }

    /* map genomic file */
    genomic_filename = gth_input_get_genomic_filename(input, gen_file_num);
    indexname = gt_str_new_cstr(genomic_filename);
    gt_str_append_char(indexname, '.');
    gt_str_append_cstr(indexname, DNASUFFIX);
    input->gen_seq_con =
      input->seq_con_constructor(gt_str_get(indexname),
                                 input->searchmode & GTHREVERSE,
                                 !translate, translate);
    gt_str_delete(indexname);
    input->genomic_translate = translate;

    /* at least one sequence in genomic virtual tree  */
    gt_assert(gth_seq_con_num_of_seqs(input->gen_seq_con) > 0);

    /* set genomic file number to new value */
    input->gen_file_num = gen_file_num;

    /* load caches, if necessary */
    if (input->use_md5_cache) {
      input->gen_md5_cache = gth_md5_cache_new(genomic_filename,
                                               input->gen_seq_con);
    }
    if (input->use_desc_cache)
      input->gen_desc_cache = gth_desc_cache_new(input->gen_seq_con);
  }
  /* else: necessary file already mapped */
  gt_assert(input->genomic_translate == translate);
}

void gth_input_load_reference_file_func(GthInput *input,
                                        unsigned long ref_file_num,
                                        bool translate,
                                        GT_UNUSED const char *src_file,
                                        GT_UNUSED int src_line)
{
  /* valid reference file number */
  gt_assert(input &&
            ref_file_num < gt_str_array_size(input->referencefiles));

#if INPUT_DEBUG
  printf("load reference (file %s, line %d): ref_file_num: %lu, translate=%s\n",
         src_file, src_line, ref_file_num, translate ? "true" : "false");
#endif

  if (input->ref_file_num != ref_file_num) {
    const char *reference_filename;
    GthAlphatype alphatype;
    GtStr *indexname;

    /* free old reference file */
    if (input->ref_file_num != GT_UNDEF_ULONG) {
      /* in this case a sequence collection has been loaded -> free it */
      gth_seq_con_delete(input->ref_seq_con);
      gth_md5_cache_delete(input->ref_md5_cache);
      /* delete caches */
      gth_desc_cache_delete(input->ref_desc_cache);
    }

    /* get alphabet type */
    alphatype = gth_input_get_alphatype(input, ref_file_num);

    /* alphabet type is valid */
    gt_assert(alphatype == DNA_ALPHA || alphatype == PROTEIN_ALPHA);

    /* loading reference sequence */
    reference_filename = gth_input_get_reference_filename(input, ref_file_num);
    indexname = gt_str_new_cstr(reference_filename);
    gt_str_append_char(indexname, '.');
    gt_str_append_cstr(indexname, alphatype == DNA_ALPHA
                                  ? DNASUFFIX
                                  : gt_str_get(input->proteinsmap));
    if (alphatype == DNA_ALPHA) {
      input->ref_seq_con = input->seq_con_constructor(gt_str_get(indexname),
                                                      true, !translate,
                                                      translate);
    }
    else {
      input->ref_seq_con = input->seq_con_constructor(gt_str_get(indexname),
                                                      false, true, true);
    }
    gt_str_delete(indexname);
    input->reference_translate = translate;

    /* at least on reference sequence in virtual tree */
    gt_assert(gth_seq_con_num_of_seqs(input->ref_seq_con) > 0);

    /* set reference file number to new value */
    input->ref_file_num = ref_file_num;

    /* load caches, if necessary */
    if (input->use_md5_cache) {
      input->ref_md5_cache = gth_md5_cache_new(reference_filename,
                                               input->ref_seq_con);
    }
    if (input->use_desc_cache)
      input->ref_desc_cache = gth_desc_cache_new(input->ref_seq_con);
  }
  /* else: necessary file already mapped */
  gt_assert(input->reference_translate == translate);
}

/* We use a ``special'' protein alphabet which has some characters which are
   wildcards in the ``normal'' protein alphabet as normal characters.
   Thereby, we can store different scores for character pairs with these
   ``wildcards''. */
static GtAlphabet* alphabet_new_blast_matrix(void)
{
  GtAlphabet *a;
  const char *blast_matrix_chars = "ARNDCQEGHILKMFPSTWYVBZX";
  char characters[2];
  size_t i, len;
  a = gt_alphabet_new_empty();
  len = strlen(blast_matrix_chars);
  characters[1] = '\0';
  for (i = 0; i < len; i++) {
    characters[0] = blast_matrix_chars[i];
    gt_alphabet_add_mapping(a, characters);
  }
  gt_alphabet_add_wildcard(a, '*');
  /* add other special wildcards */
  gt_alphabet_add_wildcard(a, 'U');
  gt_alphabet_add_wildcard(a, 'O');
  return a;
}

static GtStr* find_score_matrix_path(const char *scorematrixfile, GtError *err)
{
  GtStr *path = gt_str_new();
  int had_err = 0;
  if (gt_file_exists(scorematrixfile)) {
    gt_str_set(path, scorematrixfile);
    return path;
  }
  if (strchr(scorematrixfile, '/')) {
    gt_error_set(err, "filename \"%s\" contains illegal symbol '/': the path "
                      "list " "specified by environment variable \"%s\" cannot "
                      "be searched for it", scorematrixfile, GTHDATAENVNAME);
    had_err = -1;
  }
  if (!had_err)
    had_err = gt_file_find_in_env(path, scorematrixfile, GTHDATAENVNAME, err);
  if (!had_err && !gt_str_length(path)) {
    gt_error_set(err, "file \"%s\" not found in directory list specified by "
                 "environment variable %s", scorematrixfile, GTHDATAENVNAME);
    had_err = -1;
  }
  if (!had_err) {
    gt_assert(gt_str_length(path));
    /* path found -> append score matrix file name */
    gt_str_append_char(path, '/');
    gt_str_append_cstr(path, scorematrixfile);
  }
  if (had_err) {
    gt_str_delete(path);
    return NULL;
  }
  return path;
}

int gth_input_load_scorematrix(GthInput *input, char *scorematrixfile,
                               GthOutput *out, GtError *err)
{
  unsigned long i;
  bool protein_reffile_exists = false; /* equals true if at least one reference
                                          file has alphabet type PROTEIN_ALPHA
                                       */
  int had_err = 0;

  gt_error_check(err);

  /* loop over alphatypes to set protein_reffile_exists */
  for (i = 0; i < gt_array_size(input->alphatypes); i++) {
    if (gth_input_get_alphatype(input, i) == PROTEIN_ALPHA) {
      protein_reffile_exists = true;
      break;
    }
  }

  /* load scorematrix if necessary */
  if (protein_reffile_exists) {
    GtStr *path;

    if (out->showverbose) {
      out->showverbose("read in the following amino acid substitution matrix:");
      out->showverbose(scorematrixfile);
    }

    gt_assert(!input->score_matrix_alpha);
    input->score_matrix_alpha = alphabet_new_blast_matrix();
    if (!(path = find_score_matrix_path(scorematrixfile, err)))
      had_err = -1;
    if (!had_err) {
      input->score_matrix = gt_score_matrix_new_read(gt_str_get(path),
                                                     input->score_matrix_alpha,
                                                     err);
      if (!input->score_matrix)
        had_err = -1;
    }
    gt_str_delete(path);
  }

  return had_err;
}

GtStr* gth_input_proteinsmap(const GthInput *input)
{
  gt_assert(input);
  return input->proteinsmap;
}

GtStr* gth_input_bssmfile(const GthInput *input)
{
  gt_assert(input);
  return input->bssmfile;
}

const char* gth_input_bssmfilename(const GthInput *input)
{
  gt_assert(input);
  return gt_str_length(input->bssmfile) ? gt_str_get(input->bssmfile)
                                        : "none";
}

GtScoreMatrix* gth_input_score_matrix(const GthInput *input)
{
  gt_assert(input);
  return input->score_matrix;
}

GtAlphabet* gth_input_score_matrix_alpha(const GthInput *input)
{
  gt_assert(input);
  return input->score_matrix_alpha;
}

void gth_input_set_forward_only(GthInput *input)
{
  gt_assert(input);
  input->searchmode = GTHFORWARD;
}

void gth_input_set_reverse_only(GthInput *input)
{
  gt_assert(input);
  input->searchmode = GTHREVERSE;
}

bool gth_input_forward(const GthInput *input)
{
  gt_assert(input);
  return input->searchmode & GTHFORWARD ? true : false;
}

bool gth_input_reverse(const GthInput *input)
{
  gt_assert(input);
  return input->searchmode & GTHREVERSE ? true : false;
}

bool gth_input_both(const GthInput *input)
{
  gt_assert(input);
  return ((input->searchmode & GTHFORWARD) &&
          (input->searchmode & GTHREVERSE)) ? true : false;
}

GthAlphatype gth_input_overall_alphatype(const GthInput *input)
{
  gt_assert(input);
  return input->overall_alphatype;
}

/* the following function determines the index of the file name <ilename> in the
   array <files> and returns it. If the file name is not contained in the array,
   -1 is returned. */
static long determine_file_index(const char *filename, GtStrArray *files)
{
  unsigned long i;
  long rval = -1;
  for (i = 0; i < gt_str_array_size(files); i++) {
    if (!strcmp(filename, gt_str_array_get(files, i))) {
      rval = gt_safe_cast2long(i);
      break;
    }
  }
  return rval;
}

long gth_input_determine_genomic_file_index(const GthInput *input,
                                            const char *filename)
{
  gt_assert(input && filename);
  return determine_file_index(filename, input->genomicfiles);
}

long gth_input_determine_reference_file_index(const GthInput *input,
                                              const char *filename)
{
  gt_assert(input && filename);
  return determine_file_index(filename, input->referencefiles);
}

int gth_input_set_and_check_substring_spec(GthInput *input, GtError *err)
{
  unsigned long numofsequences, gen_total_length;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(input);

  if (input->genomicfrompos) {
    /* load genomic file first */
    /* XXX: we do not need the reverse complement here */
    gth_input_load_genomic_file(input, 0, true);

    /* now we can set numofsequences and gen_total_length */
    numofsequences   = gth_input_num_of_gen_seqs(input, 0);
    gen_total_length = gth_input_genomic_file_total_length(input, 0);

    /* at least one genomic sequence in multiseq */
    gt_assert(numofsequences > 0);

    /* position setting options have been used, check if only one genomic
       sequence is given */

    if (numofsequences > 1) {
      gt_error_set(err, "-%s only allowed for genomic file containing one "
                        "sequence", FROMPOS_OPT_CSTR);
      had_err = -1;;
    }

    /* check if the values of frompos or width are valid */
    if (!had_err && input->genomictopos > 0) {
      /* -frompos has been used together with -topos */

      /* checking if values are correctly set */
      if (input->genomicfrompos == input->genomictopos) {
        gt_error_set(err, "frompos equals topos");
        had_err = -1;
      }
      if (!had_err && input->genomicfrompos > gen_total_length) {
        gt_error_set(err, "frompos is larger than total length of genomic "
                          "sequence");
        had_err = -1;
      }
      if (!had_err && input->genomictopos > gen_total_length) {
        gt_error_set(err, "topos is larger than total length of genomic "
                          "sequence");
        had_err = -1;
      }
      if (!had_err &&
          input->genomicfrompos >= input->genomictopos) {
        gt_error_set(err, "frompos has to be less than topos");
        had_err = -1;
      }

      if (!had_err) {
        /* all values correctly set -> setting width and changing values
           because from here on we count starting with 0 */
        input->genomicfrompos -= 1;
        input->genomictopos   -= 1;
        input->genomicwidth    = input->genomictopos -
                                 input->genomicfrompos + 1;
      }
    }
    else if (!had_err) {
      /* -frompos has been used together with -width */

      /* checking if values are correctly set */
      if (input->genomicfrompos + input->genomicwidth - 1
          > gen_total_length) {
        gt_error_set(err, "frompos + width is larger than total length of "
                          "genomic sequence");
        had_err = -1;
      }

      if (!had_err) {
        /* all values correctly set -> setting topos and changing values
           because from here on we count starting with 0 */
        input->genomicfrompos -= 1;
        input->genomictopos = input->genomicfrompos +
                              input->genomicwidth - 1;
      }
    }

    input->use_substring_spec = true;
  }

  return had_err;
}

bool gth_input_use_substring_spec(const GthInput *input)
{
  gt_assert(input);
  return input->use_substring_spec;
}

unsigned long gth_input_genomic_substring_from(const GthInput *input)
{
  gt_assert(input);
  return input->genomicfrompos;
}

unsigned long gth_input_genomic_substring_to(const GthInput *input)
{
  gt_assert(input);
  return input->genomictopos;
}

unsigned long* gth_input_genomicfrompos_ptr(GthInput *input)
{
  gt_assert(input);
  return &input->genomicfrompos;
}

unsigned long* gth_input_genomicwidth_ptr(GthInput *input)
{
  gt_assert(input);
  return &input->genomicwidth;
}

unsigned long* gth_input_genomictopos_ptr(GthInput *input)
{
  gt_assert(input);
  return &input->genomictopos;
}

int gth_input_make_indices(GthInput *input, const char *progname, GtError *err)
{
  int had_err;
  GthOutput *out;
  gt_error_check(err);
  gt_assert(input);
  out = gthoutput_new();
  had_err = input->file_preprocessor(input, true, false, false, false, false,
                                     false, progname,
                                     GTH_DEFAULT_TRANSLATIONTABLE, out, err);
  if (!had_err) {
    had_err = gth_input_load_scorematrix(input, GTH_DEFAULT_SCOREMATRIX, out,
                                         err);
  }
  gthoutput_delete(out);
  return had_err;
}

void gth_input_delete_current(GthInput *input)
{
  /* free current genomic virtual tree */
  if (input->gen_file_num != GT_UNDEF_ULONG) {
    /* in this case a virtual tree has been loaded -> free it */
    gth_seq_con_delete(input->gen_seq_con);
    gth_md5_cache_delete(input->gen_md5_cache);
    gth_desc_cache_delete(input->gen_desc_cache);
  }

  /* free current reference virtual tree */
  if (input->ref_file_num != GT_UNDEF_ULONG) {
    /* in this case a virtual tree has been loaded -> free it */
    gth_seq_con_delete(input->ref_seq_con);
    gth_md5_cache_delete(input->ref_md5_cache);
    gth_desc_cache_delete(input->ref_desc_cache);
  }

  /* set the filenumbers to undefined values */
  input->gen_file_num = GT_UNDEF_ULONG;
  input->ref_file_num = GT_UNDEF_ULONG;

  input->genomic_translate = GT_UNDEF_BOOL;
  input->reference_translate = GT_UNDEF_BOOL;
}

void gth_input_delete_complete(GthInput *input)
{
  if (!input) return;
  gth_input_delete_current(input);
  gt_str_delete(input->bssmfile);
  gt_str_delete(input->proteinsmap);
  gt_score_matrix_delete(input->score_matrix);
  gt_alphabet_delete(input->score_matrix_alpha);
  gt_array_delete(input->alphatypes);
  gt_str_array_delete(input->referencefiles);
  gt_str_array_delete(input->genomicfiles);
  gt_free(input);
}

GthSeqCon* gth_input_current_gen_seq_con(GthInput *input)
{
  gt_assert(input);
  return input->gen_seq_con;
}

GthSeqCon* gth_input_current_ref_seq_con(GthInput *input)
{
  gt_assert(input);
  return input->ref_seq_con;
}

const unsigned char* gth_input_current_gen_seq_tran(const GthInput *input)
{
  gt_assert(input);
  return gth_seq_con_get_tran_seq(input->gen_seq_con, 0);
}

const unsigned char* gth_input_current_gen_seq_tran_rc(const GthInput *input)
{
  gt_assert(input);
  return gth_seq_con_get_tran_seq_rc(input->gen_seq_con, 0);
}

const unsigned char* gth_input_current_gen_seq_orig(const GthInput *input)
{
  gt_assert(input);
  return gth_seq_con_get_orig_seq(input->gen_seq_con, 0);
}

const unsigned char* gth_input_current_gen_seq_orig_rc(const GthInput *input)
{
  gt_assert(input);
  return gth_seq_con_get_orig_seq_rc(input->gen_seq_con, 0);
}

const unsigned char* gth_input_current_ref_seq_tran(const GthInput *input)
{
  gt_assert(input);
  return gth_seq_con_get_tran_seq(input->ref_seq_con, 0);
}

const unsigned char* gth_input_current_ref_seq_tran_rc(const GthInput *input)
{
  gt_assert(input);
  return gth_seq_con_get_tran_seq_rc(input->ref_seq_con, 0);
}

const unsigned char* gth_input_current_ref_seq_orig(const GthInput *input)
{
  gt_assert(input);
  return gth_seq_con_get_orig_seq(input->ref_seq_con, 0);
}

const unsigned char* gth_input_current_ref_seq_orig_rc(const GthInput *input)
{
  gt_assert(input);
  return gth_seq_con_get_orig_seq_rc(input->ref_seq_con, 0);
}

GtAlphabet* gth_input_current_gen_alphabet(GthInput *input)
{
  gt_assert(input);
  return gth_seq_con_get_alphabet(input->gen_seq_con);
}

GtAlphabet* gth_input_current_ref_alphabet(GthInput *input)
{
  gt_assert(input);
  return gth_seq_con_get_alphabet(input->ref_seq_con);
}
