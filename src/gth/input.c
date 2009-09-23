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

#include <ctype.h>
#include <math.h>
#include "core/assert_api.h"
#include "core/fileutils_api.h"
#include "core/ma.h"
#include "core/safearith.h"
#include "core/str.h"
#include "core/undef.h"
#include "gth/default.h"
#include "gth/gthdef.h"
#include "gth/input.h"
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
  GthSeqColConstructor seq_col_constructor;
  GthSeqCol *genomic_seq_col;     /* containing the preprocessed genomic
                                     sequence */
  GthSeqCol *reference_seq_col;   /* containing the reference sequence */
  GtStr *proteinsmap,             /* name of protein symbol mapping */
        *bssmfile;                /* file for bssm parameter */
  GtScoreMatrix *score_matrix;    /* contains the amino acid substitution matrix
                                     (if necessary) */
  GtAlphabet *score_matrix_alpha; /* the alphabet used for the scoring matrix */

  unsigned int searchmode;        /* stores bits GTHFORWARD and GTHREVERSE */
  bool use_substring_spec;
  unsigned long genomicfrompos,   /* analyse genomic seq. from this position */
                genomicwidth,     /* analyse only this width of genomic seq. */
                genomictopos;     /* = genomicfrompos + genomicwidth - 1 */
};

GthInput *gth_input_new(GthInputFilePreprocessor file_preprocessor,
                        GthSeqColConstructor seq_col_constructor)
{
  GthInput *input = gt_calloc(1, sizeof *input);
  gt_assert(seq_col_constructor);
  input->file_preprocessor = file_preprocessor;
  input->genomicfiles = gt_str_array_new();
  input->referencefiles = gt_str_array_new();
  input->alphatypes = gt_array_new(sizeof (GthAlphatype));
  input->overall_alphatype = UNDEF_ALPHA;
  input->gen_file_num = GT_UNDEF_ULONG;
  input->ref_file_num = GT_UNDEF_ULONG;
  input->seq_col_constructor = seq_col_constructor;
  input->proteinsmap = gt_str_new_cstr(DEFAULT_PROTEINSMAP);
  input->bssmfile = gt_str_new();
  input->searchmode = GTHFORWARD | GTHREVERSE;
  return input;
}

int gth_input_preprocess(GthInput *input,
                         bool gthconsensus,
                         bool noautoindex,
                         bool skipindexcheck,
                         bool maskpolyAtails,
                         bool online,
                         bool inverse,
                         const char *progname,
                         char *scorematrixfile,
                         unsigned int translationtable,
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

const unsigned char* gth_input_genomic_sequence(GthInput *input,
                                                unsigned long filenum,
                                                bool forward)
{
  gt_assert(input);
  gth_input_load_genomic_file(input, filenum);
  if (forward)
    return gth_seq_col_get_tran_seq(input->genomic_seq_col, 0);
  else {
    gt_assert(input->searchmode & GTHREVERSE);
    return gth_seq_col_get_tran_seq_rc(input->genomic_seq_col, 0);
  }
}

const unsigned char* gth_input_original_genomic_sequence(GthInput *input,
                                                         unsigned long filenum,
                                                         bool forward)
{
  gt_assert(input);
  gth_input_load_genomic_file(input, filenum);
  if (forward)
    return gth_seq_col_get_orig_seq(input->genomic_seq_col, 0);
  else {
    gt_assert(input->searchmode & GTHREVERSE);
    return gth_seq_col_get_orig_seq_rc(input->genomic_seq_col, 0);
  }
}

void gth_input_echo_genomic_description(GthInput *input, unsigned long filenum,
                                        unsigned long seqnum, GtFile *outfp)
{
  gt_assert(input);
  gth_input_load_genomic_file(input, filenum);
  gth_seq_col_echo_description(input->genomic_seq_col, seqnum, outfp);
}

void gth_input_echo_reference_description(GthInput *input,
                                          unsigned long filenum,
                                          unsigned long seqnum,
                                          GtFile *outfp)
{
  gt_assert(input);
  gth_input_load_reference_file(input, filenum);
  gth_seq_col_echo_description(input->reference_seq_col, seqnum, outfp);
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
                                       unsigned long filenum,
                                       unsigned long seqnum, bool forward,
                                       GtFile *outfp)
{
  unsigned char *refseq;
  unsigned long i, reflength;
  gt_assert(input);
  gth_input_load_reference_file(input, filenum);

  /* get reference sequence */
  if (forward)
    refseq = gth_seq_col_get_orig_seq(input->reference_seq_col, seqnum);
  else
    refseq = gth_seq_col_get_orig_seq_rc(input->reference_seq_col, seqnum);

  /* output reference sequence */
  reflength = gth_seq_col_get_length(input->reference_seq_col, seqnum);
  if (format)
    format_reference_seq(refseq, reflength, outfp);
  else {
    for (i = 0; i < reflength; i++)
      gt_file_xfputc(refseq[i], outfp);
  }
}

static void save_sequenceid(GtStr *sequenceid, GthSeqCol *seqcol,
                            unsigned long seqnum)
{
  unsigned long i, len;
  unsigned char *desc, cc;
  GtStr *description;

  /* sequence number is defined */
  gt_assert(seqnum != GT_UNDEF_ULONG);

  description = gt_str_new();
  gth_seq_col_get_description(seqcol, seqnum, description);
  len  = gt_str_length(description)-1; /* -1 to ignore '\n' */
  desc = (GtUchar*) gt_str_get(description);

  /* if the description ends with "\r\n" we ignore the '\r', too */
  if (len > 0 && desc[len-1] == '\r')
    len--;

  i = 0;

  if ((len >= 2) && (desc[0] == 'g') && (desc[1] == 'i') && (desc[2] == '|')) {
    /* skip 'gi|' */
    i = 3;
  }
  else if ((len >= 2) && (desc[0] == 'S') && (desc[1] == 'Q') &&
           (desc[2] == ';')) {
    /* skip 'SQ;' */
    i = 3;
  }
  else if ((len >= 3) && (desc[0] == '(') && (desc[1] == 'g') &&
           (desc[2] == 'i') && (desc[3] == '|')) {
    /* skip '(gi|' */
    i = 4;
  }
  else if ((len >= 3) && (desc[0] == 'r') && (desc[1] == 'e') &&
           (desc[2] == 'f') && (desc[3] == '|')) {
    /* skip 'ref|' */
    i = 4;
  }

  for (/* init already done */ ; i < len; i++) {
    cc = desc[i];
    if (cc == ':' || cc == '|' || cc == '\t' || cc == ' ')
      break;
    gt_str_append_char(sequenceid, desc[i]);
  }

  gt_str_delete(description);
}

void gth_input_save_gen_id(GthInput *input, GtStr *id, unsigned long file_num,
                           unsigned long seq_num)
{
  gt_assert(input && id);
  gth_input_load_genomic_file(input, file_num);
  save_sequenceid(id, input->genomic_seq_col, seq_num);
}

void gth_input_save_ref_id(GthInput *input, GtStr *id, unsigned long file_num,
                           unsigned long seq_num)
{
  gt_assert(input && id);
  gth_input_load_reference_file(input, file_num);
  save_sequenceid(id, input->reference_seq_col, seq_num);
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
                                                  unsigned long filenum)
{
  gt_assert(input);
  gth_input_load_genomic_file(input, filenum);
  return gth_seq_col_total_length(input->genomic_seq_col);
}

unsigned long gth_input_num_of_gen_seqs(GthInput *input, unsigned long filenum)
{
  gt_assert(input);
  gth_input_load_genomic_file(input, filenum);
  return gth_seq_col_num_of_seqs(input->genomic_seq_col);
}

unsigned long gth_input_num_of_ref_seqs(GthInput *input, unsigned long filenum)
{
  gt_assert(input);
  gth_input_load_reference_file(input, filenum);
  return gth_seq_col_num_of_seqs(input->reference_seq_col);
}

GtRange gth_input_get_relative_genomic_range(GthInput *input,
                                             unsigned long filenum,
                                             unsigned long seqnum)
{
  gt_assert(input);
  gth_input_load_genomic_file(input, filenum);
  return gth_seq_col_get_relative_range(input->genomic_seq_col, seqnum);
}

GtRange gth_input_get_genomic_range(GthInput *input,
                                    unsigned long filenum,
                                    unsigned long seqnum)
{
  gt_assert(input);
  gth_input_load_genomic_file(input, filenum);
  return gth_seq_col_get_range(input->genomic_seq_col, seqnum);
}

GtRange gth_input_get_reference_range(GthInput *input,
                                      unsigned long filenum,
                                      unsigned long seqnum)
{
  gt_assert(input);
  gth_input_load_reference_file(input, filenum);
  return gth_seq_col_get_range(input->reference_seq_col, seqnum);
}

void gth_input_load_genomic_file(GthInput *input, unsigned long gen_file_num)
{
  char indexname[PATH_MAX+MAXSUFFIXLEN+1];

  /* valid genomic file number */
  gt_assert(input && gen_file_num < gt_str_array_size(input->genomicfiles));

  if (input->gen_file_num != gen_file_num) {
    /* free old genomic file */
    if (input->gen_file_num != GT_UNDEF_ULONG) {
      /* in this case a sequence collection has been loaded -> free it */
      gth_seq_col_delete(input->genomic_seq_col);
    }

    /* map genomic file */
    sprintf(indexname, "%s.%s",
            gth_input_get_genomic_filename(input, gen_file_num), DNASUFFIX);
    input->genomic_seq_col =
      input->seq_col_constructor(indexname, input->searchmode & GTHREVERSE);

    /* at least one sequence in genomic virtual tree  */
    gt_assert(gth_seq_col_num_of_seqs(input->genomic_seq_col) > 0);

    /* set genomic file number to new value */
    input->gen_file_num = gen_file_num;
  }
  /* else: necessary file already mapped */
}

void gth_input_load_reference_file(GthInput *input,
                                      unsigned long ref_file_num)
{
  char indexname[PATH_MAX+MAXSUFFIXLEN+1];
  GthAlphatype alphatype;

  /* valid reference file number */
  gt_assert(input &&
            ref_file_num < gt_str_array_size(input->referencefiles));

  if (input->ref_file_num != ref_file_num) {
    /* free old reference file */
    if (input->ref_file_num != GT_UNDEF_ULONG) {
      /* in this case a sequence collection has been loaded -> free it */
      gth_seq_col_delete(input->reference_seq_col);
    }

    /* get alphabet type */
    alphatype = gth_input_get_alphatype(input, ref_file_num);

    /* alphabet type is valid */
    gt_assert(alphatype == DNA_ALPHA || alphatype == PROTEIN_ALPHA);

    /* loading reference sequence */
    sprintf(indexname, "%s.%s",
            gth_input_get_reference_filename(input, ref_file_num),
            alphatype == DNA_ALPHA ? DNASUFFIX
                                   : gt_str_get(input->proteinsmap));
    input->reference_seq_col =
      input->seq_col_constructor(indexname, alphatype == DNA_ALPHA);

    /* at least on reference sequence in virtual tree */
    gt_assert(gth_seq_col_num_of_seqs(input->reference_seq_col) > 0);

    /* set reference file number to new value */
    input->ref_file_num = ref_file_num;
  }
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

  if (!had_err) {
    /* now we can set numofsequences and gen_total_length */
    numofsequences   = gth_input_num_of_gen_seqs(input, 0);
    gen_total_length = gth_input_genomic_file_total_length(input, 0);

    /* at least one genomic sequence in multiseq */
    gt_assert(numofsequences > 0);
  }

  if (!had_err && input->genomicfrompos) {
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
                                     false, progname, DEFAULT_TRANSLATIONTABLE,
                                     out, err);
  if (!had_err)
    had_err = gth_input_load_scorematrix(input, DEFAULT_SCOREMATRIX, out, err);
  gthoutput_delete(out);
  return had_err;
}

void gth_input_delete_current(GthInput *input)
{
  /* free current genomic virtual tree */
  if (input->gen_file_num != GT_UNDEF_ULONG) {
    /* in this case a virtual tree has been loaded -> free it */
    gth_seq_col_delete(input->genomic_seq_col);
  }

  /* free current reference virtual tree */
  if (input->ref_file_num != GT_UNDEF_ULONG) {
    /* in this case a virtual tree has been loaded -> free it */
    gth_seq_col_delete(input->reference_seq_col);
  }

  /* set the filenumbers to undefined values */
  input->gen_file_num = GT_UNDEF_ULONG;
  input->ref_file_num = GT_UNDEF_ULONG;
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

GthSeqCol* gth_input_current_gen_seq_col(GthInput *input)
{
  gt_assert(input);
  return input->genomic_seq_col;
}

GthSeqCol* gth_input_current_ref_seq_col(GthInput *input)
{
  gt_assert(input);
  return input->reference_seq_col;
}

const unsigned char* gth_input_current_gen_seq_tran(const GthInput *input)
{
  gt_assert(input);
  return gth_seq_col_get_tran_seq(input->genomic_seq_col, 0);
}

const unsigned char* gth_input_current_gen_seq_tran_rc(const GthInput *input)
{
  gt_assert(input);
  return gth_seq_col_get_tran_seq_rc(input->genomic_seq_col, 0);
}

const unsigned char* gth_input_current_gen_seq_orig(const GthInput *input)
{
  gt_assert(input);
  return gth_seq_col_get_orig_seq(input->genomic_seq_col, 0);
}

const unsigned char* gth_input_current_gen_seq_orig_rc(const GthInput *input)
{
  gt_assert(input);
  return gth_seq_col_get_orig_seq_rc(input->genomic_seq_col, 0);
}

const unsigned char* gth_input_current_ref_seq_tran(const GthInput *input)
{
  gt_assert(input);
  return gth_seq_col_get_tran_seq(input->reference_seq_col, 0);
}

const unsigned char* gth_input_current_ref_seq_tran_rc(const GthInput *input)
{
  gt_assert(input);
  return gth_seq_col_get_tran_seq_rc(input->reference_seq_col, 0);
}

const unsigned char* gth_input_current_ref_seq_orig(const GthInput *input)
{
  gt_assert(input);
  return gth_seq_col_get_orig_seq(input->reference_seq_col, 0);
}

const unsigned char* gth_input_current_ref_seq_orig_rc(const GthInput *input)
{
  gt_assert(input);
  return gth_seq_col_get_orig_seq_rc(input->reference_seq_col, 0);
}

GtAlphabet* gth_input_current_gen_alphabet(GthInput *input)
{
  gt_assert(input);
  return gth_seq_col_get_alphabet(input->genomic_seq_col);
}

GtAlphabet* gth_input_current_ref_alphabet(GthInput *input)
{
  gt_assert(input);
  return gth_seq_col_get_alphabet(input->reference_seq_col);
}
