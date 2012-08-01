/*
  Copyright (c) 2004-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <expat.h>
#include "core/fa.h"
#include "core/unused_api.h"
#include "gth/intermediate.h"

#define SPLICEDALIGNMENT_TAG            "spliced_alignment"
#define REFERENCEALPHATYPE_TAG          "referencealphatype"
#define DNA_EOP_TYPE_TAG                "DNA_eop_type"
#define DNA_EOP_LENGTH_TAG              "DNA_eop_length"
#define PROTEIN_EOP_TYPE_TAG            "Protein_eop_type"
#define PROTEIN_EOP_LENGTH_TAG          "Protein_eop_length"
#define INDELCOUNT_TAG                  "indelcount"
#define GENOMICLENGTHDP_TAG             "genomiclengthDP"
#define GENOMICLENGTHTOTAL_TAG          "genomiclengthtotal"
#define GENOMICOFFSET_TAG               "genomicoffset"
#define REFERENCELENGTH_TAG             "referencelength"
#define DPSTARTPOS_TAG                  "dpstartpos"
#define DPENDPOS_TAG                    "dpendpos"
#define GENOMICFILENAME_TAG             "genomicfilename"
#define GENOMICFILEHASH_TAG             "genomicfilehash"
#define GENOMICSEQNUM_TAG               "genomicseqnum"
#define REFERENCEFILENAME_TAG           "referencefilename"
#define REFERENCEFILEHASH_TAG           "referencefilehash"
#define REFERENCESEQNUM_TAG             "referenceseqnum"
#define GENOMICID_TAG                   "genomicid"
#define REFERENCEID_TAG                 "referenceid"
#define GENOMICSTRANDISFORWARD_TAG      "genomicstrandisforward"
#define REFERENCESTRANDISFORWARD_TAG    "referencestrandisforward"
#define CUTOFFSSTART_TAG                "cutoffsstart"
#define CUTOFFSEND_TAG                  "cutoffsend"
#define GENOMICCUTOFF_TAG               "genomiccutoff"
#define REFERENCECUTOFF_TAG             "referencecutoff"
#define EOPCUTOFF_TAG                   "eopcutoff"
#define LEFTGENOMICEXONBORDER_TAG       "leftgenomicexonborder"
#define RIGHTGENOMICEXONBORDER_TAG      "rightgenomicexonborder"
#define LEFTREFERENCEEXONBORDER_TAG     "leftreferenceexonborder"
#define RIGHTREFERENCEEXONBORDER_TAG    "rightreferenceexonborder"
#define EXONSCORE_TAG                   "exonscore"
#define EXONINFO_TAG                    "exoninfo"
#define DONORSITEPROBABILITY_TAG        "donorsiteprobability"
#define ACCEPTORSITEPROBABILITY_TAG     "acceptorsiteprobability"
#define DONORSITESCORE_TAG              "donorsitescore"
#define ACCEPTORSITESCORE_TAG           "acceptorsitescore"
#define INTRONINFO_TAG                  "introninfo"
#define POLYASTART_TAG                  "polyAstart"
#define POLYAEND_TAG                    "polyAstop"
#define ALIGNMENTSCORE_TAG              "alignmentscore"
#define COVERAGE_TAG                    "coverage"
#define COVERAGEOFGENOMICSEGMENTISHIGHEST_TAG\
                                        "coverageofgenomicsegmentishighest"
#define CUMULATIVELENGTHOFSCOREDEXONS_TAG\
                                        "cumulativelengthofscoredexons"

#define ILLEGAL_DATA\
        fprintf(stderr, "illegal data in line %lu of file \"%s\"\n",\
                parseinfo->linenumber , parseinfo->outputfilename);\
        exit(EXIT_FAILURE)

#define SCANUINT\
        if (sscanf(data, "%ld", &ret) != 1 || ret < 0)\
        {\
          ILLEGAL_DATA;\
        }

#define SCANDOUBLE\
        if (sscanf(data, "%lf", &retdouble) != 1)\
        {\
          ILLEGAL_DATA;\
        }

typedef struct {
  GthSA *currentSA;
  GtStr *databuf,
        *genomicfilename,
        *referencefilename;
  unsigned long linenumber;
  const char *outputfilename;
  Eoptype eoptype;
  bool proteineop;
  Cutoffs cutoffs;
  Exoninfo exoninfo;
  Introninfo introninfo;
  GthInput *input;
  GthSAProcessFunc saprocessfunc;
  void *data;
  GtError *err;
} Parseinfo;

typedef struct {
  GthSACollection *sa_collection;
  GthSAFilter *sa_filter;
  GthStat *stat;
} SACollectionData;

static int store_in_sa_collection(void *data, GthSA *sa,
                                  GT_UNUSED const char *outputfilename,
                                  GT_UNUSED GtError *err)
{
  SACollectionData *sa_collection_data = (SACollectionData*) data;
  bool inserted;
  inserted = gth_sa_collection_insert_sa(sa_collection_data->sa_collection, sa,
                                         sa_collection_data->sa_filter,
                                         sa_collection_data->stat);
  if (!inserted) { /* unsuccessful insertion; discard sa */
    gth_sa_delete(sa);
  }
  return 0;
}

static bool parse_boolean(const char *data, Parseinfo *parseinfo)
{
  if (strcmp(data, "True") == 0)
    return true;
  else if (strcmp(data, "False") == 0)
    return false;
  ILLEGAL_DATA;
}

/*
  The following function handles the start elements of the intermediate output.
*/

static void start_element_handler(void *data, const XML_Char *name,
                                  GT_UNUSED const XML_Char **attributes)
{
  Parseinfo *parseinfo = (Parseinfo*) data;

  /* reset data buffer */
  gt_str_reset(parseinfo->databuf);

  /* perform actions depending on start tag */
  if (strcmp(name, SPLICEDALIGNMENT_TAG) == 0) {
    gt_assert(!parseinfo->currentSA);
    parseinfo->currentSA = gth_sa_new();
  }
}

/* XXX: implement this function */
static bool hashes_are_the_same(char *filehash, GT_UNUSED FILE *fp)
{
  if (!strcmp(filehash, GTH_UNDEFINED_HASH))
    return true;
  return false;
}

/* The following function processes a file. That is, it checkes if the file with
   name <filename> is already contained in the array <files>. If so, the index
   refering to this array is returned. Otherwise, the hash of the file content
   is compared with the hash <filehash>. If the hashes are the same, the
   filename is added to the array and the index is returned. Otherwise, the
   function calls exit(). */
static unsigned long process_file(GthInput *input, char *filename,
                                  char *filehash, bool isreferencefile,
                                  GthAlphatype alphatype)
{
  long fileindex;
  FILE *fp;

  if (isreferencefile)
    fileindex = gth_input_determine_reference_file_index(input, filename);
  else
    fileindex = gth_input_determine_genomic_file_index(input, filename);

  if (fileindex == -1) {
    /* file is not contained in array yet -> open file */
    fp = gt_fa_xfopen(filename, "r");

    /* check the hash */
    if (!hashes_are_the_same(filehash, fp)) {
      fprintf(stderr, "apparently file \"%s\" has changed\n", filename);
      exit(EXIT_FAILURE);
    }

    /* hashes equal -> store new file in array and return index number */
    gt_fa_xfclose(fp);
    if (isreferencefile) {
      gth_input_add_reference_file(input, filename, alphatype);
      fileindex = gth_input_num_of_ref_files(input) - 1;
    }
    else {
      gth_input_add_genomic_file(input, filename);
      fileindex = gth_input_num_of_gen_files(input) - 1;
    }
    return fileindex;
  }

  /* file is already contained in array -> return index number */
  return fileindex;
}

/*
  The following function handles the end elements of the intermediate output.
*/

static void end_element_handler(void *info, const XML_Char *name)
{
  Parseinfo *parseinfo = (Parseinfo*) info;
  GthSA *sa = parseinfo->currentSA;
  unsigned long datalength;
  double retdouble;
  long ret;
  char *data;

  /* save data and data length */
  data       = gt_str_get(parseinfo->databuf);
  datalength = gt_str_length(parseinfo->databuf);

  /* perform actions depending on end tag */
  if (strcmp(name, SPLICEDALIGNMENT_TAG) == 0) {
    /* before we store the spliced alignment we have to reverse its edit
       operations */
    gt_assert(sa && gth_sa_backtrace_path(sa));
    gth_backtrace_path_reverse(gth_sa_backtrace_path(sa));

    /* ensure that before an intron which is not in phase the edit operation
       has length 1 (only for protein spliced alignments) */
    gth_backtrace_path_ensure_length_1_before_introns(
                                                     gth_sa_backtrace_path(sa));

    if (parseinfo->saprocessfunc(parseinfo->data , sa,
                                 parseinfo->outputfilename, parseinfo->err)) {
      /* XXX */
      fprintf(stderr, "error: %s\n", gt_error_get(parseinfo->err));
      exit(EXIT_FAILURE);
    }
    /* reset current spliced alignment */
    parseinfo->currentSA = NULL;
 }
  else if (strcmp(name, REFERENCEALPHATYPE_TAG) == 0) {
    if (strcmp(data, "DNA_ALPHA") == 0)
      gth_sa_set_alphatype(sa, DNA_ALPHA);
    else if (strcmp(data, "PROTEIN_ALPHA") == 0) {
      gth_sa_set_alphatype(sa, PROTEIN_ALPHA);
    }
    else {
      ILLEGAL_DATA;
    }
  }
  else if (strcmp(name, DNA_EOP_TYPE_TAG) == 0) {
    if (strcmp(data, "match") == 0)
      parseinfo->eoptype = EOP_TYPE_MATCH;
    else if (strcmp(data, "deletion") == 0)
      parseinfo->eoptype = EOP_TYPE_DELETION;
    else if (strcmp(data, "insertion") == 0)
      parseinfo->eoptype = EOP_TYPE_INSERTION;
    else if (strcmp(data, "mismatch") == 0)
      parseinfo->eoptype = EOP_TYPE_MISMATCH;
    else if (strcmp(data, "intron") == 0)
      parseinfo->eoptype = EOP_TYPE_INTRON;
    else {
      ILLEGAL_DATA;
    }
  }
  else if (strcmp(name, DNA_EOP_LENGTH_TAG) == 0) {
    SCANUINT;
    gth_backtrace_path_add_eop(gth_sa_backtrace_path(sa), parseinfo->eoptype,
                               ret);
  }
  else if (strcmp(name, PROTEIN_EOP_TYPE_TAG) == 0) {
    if (strcmp(data, "match") == 0)
      parseinfo->eoptype = EOP_TYPE_MATCH;
    else if (strcmp(data, "deletion") == 0)
      parseinfo->eoptype = EOP_TYPE_DELETION;
    else if (strcmp(data, "insertion") == 0)
      parseinfo->eoptype = EOP_TYPE_INSERTION;
    else if (strcmp(data, "mismatch") == 0)
      parseinfo->eoptype = EOP_TYPE_MISMATCH;
    else if (strcmp(data, "intron") == 0)
      parseinfo->eoptype = EOP_TYPE_INTRON;
    else if (strcmp(data, "mismatch_with_1_gap") == 0)
      parseinfo->eoptype = EOP_TYPE_MISMATCH_WITH_1_GAP;
    else if (strcmp(data, "mismatch_with_2_gaps") == 0)
      parseinfo->eoptype = EOP_TYPE_MISMATCH_WITH_2_GAPS;
    else if (strcmp(data, "deletion_with_1_gap") == 0)
      parseinfo->eoptype = EOP_TYPE_DELETION_WITH_1_GAP;
    else if (strcmp(data, "deletion_with_2_gaps") == 0)
      parseinfo->eoptype = EOP_TYPE_DELETION_WITH_2_GAPS;
    else if (strcmp(data, "intron_with_1_base_left") == 0)
      parseinfo->eoptype = EOP_TYPE_INTRON_WITH_1_BASE_LEFT;
    else if (strcmp(data, "intron_with_2_bases_left") == 0)
      parseinfo->eoptype = EOP_TYPE_INTRON_WITH_2_BASES_LEFT;
    else {
      ILLEGAL_DATA;
    }
  }
  else if (strcmp(name, PROTEIN_EOP_LENGTH_TAG) == 0) {
    SCANUINT;
    gth_backtrace_path_add_eop(gth_sa_backtrace_path(sa), parseinfo->eoptype,
                               ret);
  }
  else if (strcmp(name, INDELCOUNT_TAG) == 0) {
    SCANUINT;
    /* ignore indelcount, gets recomputed anyway */
  }
  else if (strcmp(name, GENOMICLENGTHDP_TAG) == 0) {
    SCANUINT;
    gth_sa_set_gen_dp_length(sa, ret);
  }
  else if (strcmp(name, GENOMICLENGTHTOTAL_TAG) == 0) {
    SCANUINT;
    gth_sa_set_gen_total_length(sa, ret);
  }
  else if (strcmp(name, GENOMICOFFSET_TAG) == 0) {
    SCANUINT;
    gth_sa_set_gen_offset(sa, ret);
  }
  else if (strcmp(name, REFERENCELENGTH_TAG) == 0) {
    SCANUINT;
    gth_sa_set_ref_total_length(sa, ret);
  }
  else if (strcmp(name, DPSTARTPOS_TAG) == 0) {
    SCANUINT;
    gth_sa_set_gen_dp_start(sa, ret);
  }
  else if (strcmp(name, DPENDPOS_TAG) == 0) {
    SCANUINT;
    /* ignore DP end pos, gets recomputed from gen_dp_length anyway */
    gt_assert(gth_sa_gen_dp_end(sa) == ret);
  }
  else if (strcmp(name, GENOMICFILENAME_TAG) == 0) {
    /* save genomic file name */
    gt_str_append_cstr_nt(parseinfo->genomicfilename, data, datalength);
  }
  else if (strcmp(name, GENOMICFILEHASH_TAG) == 0) {
    gth_sa_set_gen_file_num(sa, process_file(parseinfo->input,
                            gt_str_get(parseinfo->genomicfilename), data, false,
                            UNDEF_ALPHA));
    /* reset genomic filename */
    gt_str_reset(parseinfo->genomicfilename);
  }
  else if (strcmp(name, GENOMICSEQNUM_TAG) == 0) {
    SCANUINT;
    gth_sa_set_gen_seq_num(sa, ret);
  }
  else if (strcmp(name, REFERENCEFILENAME_TAG) == 0) {
    /* save reference file name */
    gt_str_append_cstr_nt(parseinfo->referencefilename, data, datalength);
  }
  else if (strcmp(name, REFERENCEFILEHASH_TAG) == 0) {
    gth_sa_set_ref_file_num(sa, process_file(parseinfo->input,
                                       gt_str_get(parseinfo->referencefilename),
                                                  data, true,
                                                  gth_sa_alphatype(sa)));

    /* reset reference filename */
    gt_str_reset(parseinfo->referencefilename);
  }
  else if (strcmp(name, REFERENCESEQNUM_TAG) == 0) {
    SCANUINT;
    gth_sa_set_ref_seq_num(sa, ret);
  }
  else if (strcmp(name, GENOMICID_TAG) == 0)
    gth_sa_set_gen_id(sa, data);
  else if (strcmp(name, REFERENCEID_TAG) == 0)
    gth_sa_set_ref_id(sa, data);
  else if (strcmp(name, GENOMICSTRANDISFORWARD_TAG) == 0)
    gth_sa_set_gen_strand(sa, parse_boolean(data, parseinfo));
  else if (strcmp(name, REFERENCESTRANDISFORWARD_TAG) == 0)
    gth_sa_set_ref_strand(sa, parse_boolean(data, parseinfo));
  else if (strcmp(name, GENOMICCUTOFF_TAG) == 0) {
    SCANUINT;
    parseinfo->cutoffs.genomiccutoff = ret;
  }
  else if (strcmp(name, REFERENCECUTOFF_TAG) == 0) {
    SCANUINT;
    parseinfo->cutoffs.referencecutoff = ret;
  }
  else if (strcmp(name, EOPCUTOFF_TAG) == 0) {
    SCANUINT;
    parseinfo->cutoffs.eopcutoff = ret;
  }
  else if (strcmp(name, CUTOFFSSTART_TAG) == 0)
    gth_sa_set_cutoffs_start(sa, &parseinfo->cutoffs);
  else if (strcmp(name, CUTOFFSEND_TAG) == 0)
    gth_sa_set_cutoffs_end(sa, &parseinfo->cutoffs);
  else if (strcmp(name, LEFTGENOMICEXONBORDER_TAG) == 0) {
    SCANUINT;
    parseinfo->exoninfo.leftgenomicexonborder = ret;
  }
  else if (strcmp(name, RIGHTGENOMICEXONBORDER_TAG) == 0) {
    SCANUINT;
    parseinfo->exoninfo.rightgenomicexonborder = ret;
  }
  else if (strcmp(name, LEFTREFERENCEEXONBORDER_TAG) == 0) {
    SCANUINT;
    parseinfo->exoninfo.leftreferenceexonborder = ret;
  }
  else if (strcmp(name, RIGHTREFERENCEEXONBORDER_TAG) == 0) {
    SCANUINT;
    parseinfo->exoninfo.rightreferenceexonborder = ret;
  }
  else if (strcmp(name, EXONSCORE_TAG) == 0) {
    SCANDOUBLE;
    parseinfo->exoninfo.exonscore = retdouble;
  }
  else if (strcmp(name, EXONINFO_TAG) == 0)
    gth_sa_add_exon(sa, &parseinfo->exoninfo);
  else if (strcmp(name, DONORSITEPROBABILITY_TAG) == 0) {
    SCANDOUBLE;
    parseinfo->introninfo.donorsiteprobability = (GthFlt) retdouble;
  }
  else if (strcmp(name, ACCEPTORSITEPROBABILITY_TAG) == 0) {
    SCANDOUBLE;
    parseinfo->introninfo.acceptorsiteprobability = (GthFlt) retdouble;
  }
  else if (strcmp(name, DONORSITESCORE_TAG) == 0) {
    SCANDOUBLE;
    parseinfo->introninfo.donorsitescore = retdouble;
  }
  else if (strcmp(name, ACCEPTORSITESCORE_TAG) == 0) {
    SCANDOUBLE;
    parseinfo->introninfo.acceptorsitescore = retdouble;
  }
  else if (strcmp(name, INTRONINFO_TAG) == 0)
    gth_sa_add_intron(sa, &parseinfo->introninfo);
  else if (strcmp(name, POLYASTART_TAG) == 0) {
    SCANUINT;
    gth_sa_set_polyAtail_start(sa, ret);
  }
  else if (strcmp(name, POLYAEND_TAG) == 0) {
    SCANUINT;
    gth_sa_set_polyAtail_stop(sa, ret);
  }
  else if (strcmp(name, ALIGNMENTSCORE_TAG) == 0) {
    SCANDOUBLE;
    gth_sa_set_score(sa, retdouble);
  }
  else if (strcmp(name, COVERAGE_TAG) == 0) {
    SCANDOUBLE;
    gth_sa_set_coverage(sa, retdouble);
  }
  else if (strcmp(name, COVERAGEOFGENOMICSEGMENTISHIGHEST_TAG) == 0) {
    gth_sa_set_highest_cov(sa, parse_boolean(data, parseinfo));
  }
  else if (strcmp(name, CUMULATIVELENGTHOFSCOREDEXONS_TAG) == 0) {
    SCANUINT;
    gth_sa_set_cumlen_scored_exons(sa, ret);
  }
}

/*
  The following function handles the end charachter data of the intermediate
  output.
*/

static void character_data_handler(void *data, const XML_Char *string, int len)
{
  Parseinfo *parseinfo = (Parseinfo*) data;
  /* add data to the data buffer */
  gt_str_append_cstr_nt(parseinfo->databuf, string, len);
}

int gt_parse_intermediate_output(GthInput *input,
                                 GthSAProcessFunc saprocessfunc, void *data,
                                 const char *outputfilename,
                                 GtFile *intermediate_fp, GtError *err)
{
  GtStr *line;
  XML_Parser parser;
  enum XML_Error error;
  Parseinfo parseinfo;
  int had_err = 0;

  gt_error_check(err);

  /* init */
  line = gt_str_new();

  /* create parser */
  parser = XML_ParserCreate(NULL);

  /* init parse info structure */
  parseinfo.currentSA         = NULL;
  parseinfo.databuf           = gt_str_new();
  parseinfo.genomicfilename   = gt_str_new();
  parseinfo.referencefilename = gt_str_new();
  parseinfo.linenumber        = 0;
  parseinfo.outputfilename    = outputfilename;
  parseinfo.input             = input;
  parseinfo.saprocessfunc     = saprocessfunc;
  parseinfo.data              = data;
  parseinfo.err               = err;

  /* set element handler */
  XML_SetElementHandler(parser, start_element_handler, end_element_handler);

  /* set character data handler */
  XML_SetCharacterDataHandler(parser, character_data_handler);

  /* register the parse info structure */
  XML_SetUserData(parser, &parseinfo);

  /* parse the intermediate output line by line */
  while (gt_str_read_next_line_generic(line, intermediate_fp) != EOF) {
    parseinfo.linenumber++;

    if (XML_Parse(parser, gt_str_get(line), gt_str_length(line), false) ==
        XML_STATUS_ERROR) {
      error = XML_GetErrorCode(parser);
      gt_error_set(err, "an error occured parsing line %lu of file \"%s\": %s",
                   parseinfo.linenumber, outputfilename,
                   XML_ErrorString(error));
      had_err = -1;
    }

    /* reset line buffer */
    gt_str_reset(line);
  }

  if (!had_err) {
    /* finish parsing */
    if (XML_Parse(parser, NULL, 0, true) == XML_STATUS_ERROR) {
      error = XML_GetErrorCode(parser);
      gt_error_set(err, "an error occured while finishing the parsing of file "
                        "\"%s\": %s", outputfilename, XML_ErrorString(error));
      had_err = -1;
    }
  }

  /* free space */
  XML_ParserFree(parser);
  gt_str_delete(line);
  gt_str_delete(parseinfo.databuf);
  gt_str_delete(parseinfo.genomicfilename);
  gt_str_delete(parseinfo.referencefilename);

  return had_err;
}

bool gth_intermediate_output_is_correct(char *outputfilename,
                                        GthSACollection *orig_sa_collection,
                                        GthInput *input,
                                        GtFile **outfp, GtError *err)
{
  SACollectionData sa_collection_data;
  GthSACollection *read_sa_collection;
  GtFileMode file_mode;
  bool rval;
#ifndef NDEBUG
  unsigned long numofgenomicfiles, numofreferencefiles;
#endif

  gt_error_check(err);
  gt_assert(outputfilename);
  gt_assert(*outfp);
#ifndef NDEBUG
  numofgenomicfiles   = gth_input_num_of_gen_files(input);
  numofreferencefiles = gth_input_num_of_ref_files(input);
#endif

  /* init */
  read_sa_collection = gth_sa_collection_new(GTH_DC_NONE);
  sa_collection_data.sa_collection = read_sa_collection;
  sa_collection_data.sa_filter = NULL;
  sa_collection_data.stat = NULL;

  /* store file mode */
  file_mode = gt_file_mode(*outfp);

  /* close output file */
  gt_file_delete(*outfp);

  /* open intermediate file again for reading */
  *outfp = gt_file_xopen_file_mode(file_mode, outputfilename, "r");
  gt_assert(*outfp);

  /* read in the intermediate output */
  if (gt_parse_intermediate_output(input, store_in_sa_collection,
                                   &sa_collection_data, outputfilename, *outfp,
                                   err)) {
    fprintf(stderr, "error: %s\n", gt_error_get(err));
    exit(EXIT_FAILURE);
  }

  /* array of genomic files did not grow */
  gt_assert(numofgenomicfiles == gth_input_num_of_gen_files(input));
  /* array of reference files did not grow */
  gt_assert(numofreferencefiles == gth_input_num_of_ref_files(input));

  /* compare the trees */
  rval = gth_sa_collections_are_equal(orig_sa_collection, read_sa_collection);

  /* free */
  gth_sa_collection_delete(read_sa_collection);

  return rval;
}

static void show_parse_file_status(GthShowVerbose showverbose,
                                   unsigned long filenum,
                                   unsigned long numoffiles,
                                   const char *filename)
{
  GtStr *buf = gt_str_new();
  gt_str_append_cstr(buf, "process file ");
  gt_str_append_ulong(buf, filenum + 1);
  gt_str_append_char(buf, '/');
  gt_str_append_ulong(buf, numoffiles);
  gt_str_append_cstr(buf, ": ");
  gt_str_append_cstr(buf, filename);
  showverbose(gt_str_get(buf));
  gt_str_delete(buf);
}

/* The following function processes a set of consensus files. If no consensus
   file is given, stdin is used as input. */

int gth_process_intermediate_files(GthInput *input, GtStrArray *consensusfiles,
                                   GthSAProcessFunc saprocessfunc, void *data,
                                   GthShowVerbose showverbose, GtError *err)
{
  unsigned long i;
  GtFile *fp, *genfile;
  int had_err = 0;

  gt_error_check(err);

  /* process all files */
  if (gt_str_array_size(consensusfiles)) {
    for (i = 0; !had_err && i < gt_str_array_size(consensusfiles); i++) {
      /* open file */
      fp = gt_file_xopen(gt_str_array_get(consensusfiles, i), "r");

      if (showverbose) {
        show_parse_file_status(showverbose, i,
                               gt_str_array_size(consensusfiles),
                               gt_str_array_get(consensusfiles, i));
      }

      had_err = gt_parse_intermediate_output(input, saprocessfunc, data,
                                          gt_str_array_get(consensusfiles, i),
                                          fp, err);

      /* close file */
      gt_file_delete(fp);
    }
  }
  else {
    genfile = gt_file_new_from_fileptr(stdin);
    had_err = gt_parse_intermediate_output(input, saprocessfunc, data, "stdin",
                                           genfile, err);
    gt_file_delete_without_handle(genfile);
  }

  return had_err;
}

int gth_build_sa_collection(GthSACollection *sa_collection, GthInput *input,
                            GtStrArray *consensusfiles, GthSAFilter *sa_filter,
                            GthStat *stat, GthShowVerbose showverbose,
                            GtError *err)
{
  SACollectionData sa_collection_data;
  gt_error_check(err);
  sa_collection_data.sa_collection = sa_collection;
  sa_collection_data.sa_filter = sa_filter;
  sa_collection_data.stat = stat;

  return gth_process_intermediate_files(input, consensusfiles,
                                        store_in_sa_collection,
                                        &sa_collection_data, showverbose, err);
}
