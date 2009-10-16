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

#include "core/assert_api.h"
#include "core/cstr.h"
#include "core/error.h"
#include "core/mailaddress.h"
#include "core/outputfile.h"
#include "core/undef.h"
#include "core/unused_api.h"
#include "gth/default.h"
#include "gth/gthdef.h"
#include "gth/call_info.h"
#include "gth/gthalphatype.h"
#include "gth/gthspeciestab.h"
#include "gth/parse_options.h"
#include "gth/stat.h"

#define SHOWINTRONMAXLEN_OPT_CSTR   "showintronmaxlen"
#define TOPOS_OPT_CSTR              "topos"
#define WIDTH_OPT_CSTR              "width"
#define ICINITIALDELTA_OPT_CSTR     "icinitialdelta"
#define DPMININTRONLENGTH_OPT_CSTR  "dpminintronlen"
#define MINMATCHLEN_OPT_CSTR        "minmatchlen"
#define SEEDLENGTH_OPT_CSTR         "seedlength"
#define PRMINMATCHLEN_OPT_CSTR      "prminmatchlen"
#define PRSEEDLENGTH_OPT_CSTR       "prseedlength"
#define PROTEINSMAP_OPT_CSTR        "proteinsmap"

#define VM_MAXNUMBEROFFILES         256

/* The following function determines the species number for a given
   species name. If the species name does not occur in the table, then an
   errror is returned. */
static int findspeciesnum(char *specname, GtError *err)
{
  unsigned int i;
  gt_error_check(err);
  for (i = 0; i < NUMOFSPECIES; i++) {
    if (!strcmp(speciestab[i],specname))
      return i;
  }
  gt_error_set(err, "illegal speciesname \"%s\"", specname);
  return -1;
}

static const char *cutoff_modes[] = { "RELAXED", "STRICT", "MINIMAL" };

static int get_cutoffs_mode_from_table(char *searchstring)
{
  unsigned int i;
  for (i = 0; i < NUMOFCUTOFFMODES; i++) {
    if (!strcmp(searchstring, cutoff_modes[i]))
      return  i;
  }
  gt_assert(0); /* cannot happen */
  return NUMOFCUTOFFMODES; /* shut up compiler */
}

static int show_gth_help_trailer(GT_UNUSED const char *progname,
                                 GT_UNUSED void *data, GT_UNUSED GtError *err)
{
  gt_error_check(err);
  printf("\nFor detailed information, please refer to the manual of "
         "GenomeThreader.");
  return 0;
}

GtOPrval gth_parse_options(GthCallInfo *callinfo, GthInput *input,
                           int *parsed_args, int argc, const char **argv,
                           bool gthconsensus_parsing,
                           GtStrArray *consensusfiles, GthStat *stat,
                           void(*showverbose)(const char *),
                           void(*showverboseVM)(char *),
                           GtShowVersionFunc show_version, GtError *err)
{
  unsigned long i;
  int ret, mode;
  double u12donorprob, u12donorprob1mism;
  GtStrArray *genomic_files, *cdna_files, *protein_files;
  GtStr *specbuf, *parsed_species, *chaining_param, *leadcutoffsmode,
        *termcutoffsmode;
  bool forward = false, reverse = false, verbose, mincutoffs = false,
       nou12intronmodel, exondistri, introndistri, refseqcovdistri = false,
       matchnumdistri = false;
  GtOptionParser *op;
  GtOutputFileInfo *ofi;
  GtOption *optgenomic = NULL,              /* mandatory input */
         *optcdna = NULL,                 /* mandatory input */
         *optprotein = NULL,              /* mandatory input */
         *optspecies = NULL,              /* BSSM parameter */
         *optbssm = NULL,                 /* BSSM parameter */
         *optscorematrix = NULL,          /* score matrix */
         *opttranslationtable = NULL,     /* translationtable */
         *optforward = NULL,              /* strand direction */
         *optreverse = NULL,              /* strand direction */
         *optfrompos = NULL,              /* genomic sequence positions */
         *opttopos = NULL,                /* genomic sequence positions */
         *optwidth = NULL,                /* genomic sequence positions */
         *optverbose = NULL,              /* output */
         *optverboseseqs = NULL,          /* output */
         *optcomments = NULL,             /* output */
         *optxmlout = NULL,               /* output */
         *optgff3out = NULL,              /* output */
         *optskipalignmentout = NULL,     /* output */
         *optmincutoffs = NULL,           /* output */
         *optshowintronmaxlen = NULL,     /* output */
         *optminorflength = NULL,         /* output */
         *optshowseqnums = NULL,          /* output */
         *optgs2out = NULL,               /* output */
         *optmaskpolyatails = NULL,       /* data preprocessing */
         *optproteinsmap = NULL,          /* data preprocessing */
         *optnoautoindex = NULL,          /* data preprocessing */
         *optcreateindicesonly = NULL,    /* data preprocessing */
         *optskipindexcheck = NULL,       /* data preprocessing */
         *optminmatchlen = NULL,          /* sim. filter, vmatch, dna matching
                                           */
         *optseedlength = NULL,           /* sim. filter, vmatch, dna matching
                                           */
         *optexdrop = NULL,               /* sim. filter, vmatch, dna matching
                                           */
         *optprminmatchlen = NULL,        /* sim. filter, vmatch, pr. matching
                                           */
         *optprseedlength = NULL,         /* sim. filter, vmatch, pr. matching
                                           */
         *optprhdist = NULL,              /* sim. filter, vmatch, pr. matching
                                           */
         *optonline = NULL,               /* sim. filter, vmatch */
         *optinverse = NULL,              /* sim. filter, vmatch */
         *optexact = NULL,                /* sim. filter, vmatch */
         *optedist = NULL,                /* sim. filter, vmatch */
         *optmaxnumofmatches = NULL,      /* sim. filter, vmatch */
         *optfragweightfactor = NULL,     /* sim. filter, before gl. chaining */
         *optgcmaxgapwidth = NULL,        /* sim. filter, global chaining */
         *optrare = NULL,                 /* sim. filter, global chaining */
         *optparalogs = NULL,             /* sim. filter, global chaining */
         *optenrichchains = NULL,         /* sim. filter, global chaining */
         *optgcmincoverage = NULL,        /* sim. filter, gl. chaining filter */
         *optstopafterchaining = NULL,    /* sim. filter, gl. chaining filter */
         *optintroncutout = NULL,         /* sim. filter, after gl. chaining */
         *optautointroncutout = NULL,     /* sim. filter, after gl. chaining */
         *opticinitialdelta = NULL,       /* sim. filter, after gl. chaining */
         *opticiterations = NULL,         /* sim. filter, after gl. chaining */
         *opticdeltaincrease = NULL,      /* sim. filter, after gl. chaining */
         *opticminremlength = NULL,       /* sim. filter, after gl. chaining */
         *optnou12intronmodel = NULL,     /* U12-type intron model */
         *optu12donorprob = NULL,         /* U12-type intron model */
         *optu12donorprob1mism = NULL,    /* U12-type intron model */
         *optprobies = NULL,              /* basic DP algorithm */
         *optprobdelgen = NULL,           /* basic DP algorithm */
         *optidentityweight = NULL,       /* basic DP algorithm */
         *optmismatchweight = NULL,       /* basic DP algorithm */
         *optundetcharweight = NULL,      /* basic DP algorithm */
         *optdeletionweight = NULL,       /* basic DP algorithm */
         *optfreeintrontrans = NULL,      /* basic DP algorithm */
         *optdpminexonlength = NULL,      /* short exon/intron parameters */
         *optdpminintronlength = NULL,    /* short exon/intron parameters */
         *optshortexonpenalty = NULL,     /* short exon/intron parameters */
         *optshortintronpenalty = NULL,   /* short exon/intron parameters */
         *optwzerotransition = NULL,      /* special parameters for DP
                                             algorithm */
         *optwdecreasedoutput = NULL,     /* special parameters for DP
                                             algorithm */
         *optdetectsmallexons = NULL,     /* DP (detect small exons) */
         *optnoicinintroncheck = NULL,    /* DP (intron cutout technique) */
         *optproteinexonpenal = NULL,     /* protein DP */
         *optleadcutoffsmode = NULL,      /* proc. of `raw' spliced
                                             alignments */
         *opttermcutoffsmode = NULL,      /* proc. of `raw' spliced
                                             alignments */
         *optcutoffsminexonlength = NULL, /* proc. of `raw' spliced
                                             alignments */
         *optscoreminexonlength = NULL,   /* proc. of `raw' spliced
                                             alignments */
         *optminaveragessp = NULL,        /* advanced similarity filter option
                                           */
         *optintermediate = NULL,         /* stop after SA computation */
         *optsortags = NULL,              /* postproc. of PGLs, sorting of
                                             AGSs */
         *optsortagswf = NULL,            /* postproc. of PGLs, sorting of
                                             AGSs */
         *optdisableclustersas = NULL,    /* consensus phase
                                             (development option) */
         *optcdnaforwardonly = NULL,      /* similarity filter
                                             (development option) */
         *optexondistri = NULL,           /* statistics */
         *optintrondistri = NULL,         /* statistics */
         *optrefseqcovdistri = NULL,      /* statistics */
         *optmatchnumdistri = NULL,       /* statistics */
         *optfirstalshown = NULL,         /* miscellaneous */
         *optshoweops = NULL;             /* testing */
  GtOPrval oprval;

  gt_error_check(err);
  gt_assert(callinfo && input);

  /* init */
  if (gthconsensus_parsing) {
    op = gt_option_parser_new("[option ...] [file ...]", "Show GenomeThreader "
                              "output files containing intermediate results "
                              "and assemble\nthe contained spliced alignments "
                              "to consensus spliced alignments.");
  }
  else {
    op = gt_option_parser_new("[option ...] -genomic file [...] -cdna file "
                              "[...] -protein file [...]",
                              "Compute similarity-based gene structure "
                              "predictions (spliced alignments)\n"
                              "using cDNA/EST and/or protein sequences and "
                              "assemble the resulting spliced\n"
                              "alignments to consensus spliced alignments.");
  }

  ofi = gt_outputfileinfo_new();
  genomic_files   = gt_str_array_new();
  cdna_files      = gt_str_array_new();
  protein_files   = gt_str_array_new();
  chaining_param  = gt_str_new();
  leadcutoffsmode = gt_str_new();
  termcutoffsmode = gt_str_new();

  /* -genomic */
  if (!gthconsensus_parsing) {
    optgenomic = gt_option_new_filenamearray("genomic", "specify input files "
                                          "containing genomic sequences\n"
                                          "mandatory option", genomic_files);
    gt_option_parser_add_option(op, optgenomic);
  }

  /* -cdna */
  if (!gthconsensus_parsing) {
    optcdna = gt_option_new_filenamearray("cdna", "specifiy input files "
                                       "containing cDNA/EST sequences",
                                       cdna_files);
    gt_option_parser_add_option(op, optcdna);
  }

  /* -protein */
  if (!gthconsensus_parsing) {
    optprotein = gt_option_new_filenamearray("protein", "specify input files "
                                          "containing protein sequences",
                                          protein_files);
    gt_option_parser_add_option(op, optprotein);
  }

  /* fill <specbuf> with the names of the species */
  specbuf = gt_str_new_cstr("specify species to select splice site model which "
                            "is most appropriate; possible species:\n");
  for (i = 0; i < NUMOFSPECIES; i++) {
    gt_str_append_char(specbuf, '"');
    gt_str_append_cstr(specbuf, speciestab[i]);
    gt_str_append_char(specbuf, '"');
    if (i < NUMOFSPECIES - 1)
      gt_str_append_char(specbuf, '\n');
  }

  /* -species */
  parsed_species = gt_str_new();
  if (!gthconsensus_parsing) {
    optspecies = gt_option_new_string("species", gt_str_get(specbuf),
                                      parsed_species, NULL);

    gt_option_parser_add_option(op, optspecies);
  }

  /* -bssm */
  if (!gthconsensus_parsing) {
    optbssm = gt_option_new_string("bssm", "read bssm parameter from file in "
                                   "the path given by the environment variable "
                                   "BSSMDIR", gth_input_bssmfile(input), NULL);
    gt_option_parser_add_option(op, optbssm);
  }

  /* -scorematrix */
  optscorematrix = gt_option_new_string("scorematrix", "read amino acid "
                                     "substitution scoring matrix from file in "
                                     "the path given by the environment "
                                      "variable GTHDATADIR",
                                      callinfo->scorematrixfile,
                                      DEFAULT_SCOREMATRIX);
  if (gthconsensus_parsing)
    gt_option_is_extended_option(optscorematrix);
  gt_option_parser_add_option(op, optscorematrix);

  /* -translationtable */
  opttranslationtable = gt_option_new_uint("translationtable", "set the codon "
                                           "translation table used for codon "
                                           "translation in matching, DP, and "
                                           "output",
                                           &callinfo->translationtable,
                                           DEFAULT_TRANSLATIONTABLE);
  gt_option_parser_add_option(op, opttranslationtable);

  /* -f */
  if (!gthconsensus_parsing) {
    optforward = gt_option_new_bool("f", "analyze only forward strand of "
                                    "genomic sequences", &forward, false);
    gt_option_parser_add_option(op, optforward);
  }

  /* -r */
  if (!gthconsensus_parsing) {
    optreverse = gt_option_new_bool("r", "analyze only reverse strand of "
                                    "genomic sequences", &reverse, false);
    gt_option_parser_add_option(op, optreverse);
  }

  /* -frompos */
  if (!gthconsensus_parsing) {
    optfrompos = gt_option_new_ulong_min(FROMPOS_OPT_CSTR, "analyze genomic "
                                      "sequence from this position\nrequires "
                                      "-topos or -width; counting from 1 on",
                                     gth_input_genomicfrompos_ptr(input),
                                      0, 1);
    gt_option_parser_add_option(op, optfrompos);
  }

  /* -topos */
  if (!gthconsensus_parsing) {
    opttopos = gt_option_new_ulong_min(TOPOS_OPT_CSTR, "analyze genomic "
                                       "sequence to this position\nrequires "
                                       "-frompos; counting from 1 on",
                                       gth_input_genomictopos_ptr(input), 0, 1);
    gt_option_parser_add_option(op, opttopos);
  }

  /* -width */
  if (!gthconsensus_parsing) {
    optwidth = gt_option_new_ulong_min(WIDTH_OPT_CSTR, "analyze only this "
                                       "width of genomic sequence\nrequires "
                                       "-frompos",
                                       gth_input_genomicwidth_ptr(input), 0, 1);
    gt_option_parser_add_option(op, optwidth);
  }

  /* -v */
  optverbose = gt_option_new_verbose(&verbose);
  gt_option_parser_add_option(op, optverbose);

  /* -verboseseqs */
  optverboseseqs = gt_option_new_bool("verboseseqs", "show additonal sequence "
                                   "information\n(for debugging purposes)",
                                   &callinfo->out->verboseseqs, false);
  gt_option_is_development_option(optverboseseqs);
  gt_option_parser_add_option(op, optverboseseqs);

  /* -comments */
  optcomments = gt_option_new_bool("comments", "show (additional) comments",
                                &callinfo->out->comments, false);
  gt_option_is_extended_option(optcomments);
  gt_option_is_development_option(optcomments);
  gt_option_parser_add_option(op, optcomments);

  /* -xmlout */
  optxmlout = gt_option_new_bool("xmlout", "show output in XML format",
                              &callinfo->out->xmlout, false);
  gt_option_parser_add_option(op, optxmlout);

  /* -gff3out */
  optgff3out = gt_option_new_bool("gff3out" , "show output in GFF3 format",
                               &callinfo->out->gff3out, false);
  gt_option_is_development_option(optgff3out);
  gt_option_parser_add_option(op, optgff3out);

  /* output file options */
  gt_outputfile_register_options(op, &callinfo->out->outfp, ofi);

  /* -skipalignmentout */
  optskipalignmentout = gt_option_new_bool("skipalignmentout", "skip output of "
                                        "spliced alignments",
                                        &callinfo->out->skipalignmentout,
                                        false);
  gt_option_is_extended_option(optskipalignmentout);
  gt_option_parser_add_option(op, optskipalignmentout);

  /* -mincutoffs */
  if (!gthconsensus_parsing) {
    optmincutoffs = gt_option_new_bool("mincutoffs", "show full spliced "
                                    "alignments\ni.e., cutoffs mode for "
                                    "leading and terminal bases is MINIMAL",
                                    &mincutoffs, false);
    gt_option_is_extended_option(optmincutoffs);
    gt_option_parser_add_option(op, optmincutoffs);
  }

  /* -showintronmaxlen */
  optshowintronmaxlen = gt_option_new_ulong(SHOWINTRONMAXLEN_OPT_CSTR, "set "
                                            "the maximum length of a fully "
                                            "shown intron\nIf set to 0, all "
                                            "introns are shown completely",
                                            &callinfo->out->showintronmaxlen,
                                            DEFAULT_SHOWINTRONMAXLEN);
  gt_option_is_extended_option(optshowintronmaxlen);
  gt_option_parser_add_option(op, optshowintronmaxlen);

  /* -minorflen */
  optminorflength = gt_option_new_ulong_min("minorflen", "set the minimum "
                                            "length of an ORF to be shown",
                                            &callinfo->out->minORFlength,
                                            DEFAULT_MINORFLENGTH, 1);
  gt_option_is_extended_option(optminorflength);
  gt_option_parser_add_option(op, optminorflength);

  /* -showseqnums */
  optshowseqnums = gt_option_new_bool("showseqnums", "show sequence numbers in "
                                   "output", &callinfo->out->showseqnums,
                                   DEFAULT_SHOWSEQNUMS);
  gt_option_is_extended_option(optshowseqnums);
  gt_option_parser_add_option(op, optshowseqnums);

  /* -gs2out */
  optgs2out = gt_option_new_bool("gs2out", "output in old GeneSeqer2 format",
                              &callinfo->out->gs2out, DEFAULT_GS2OUT);
  gt_option_parser_add_option(op, optgs2out);

  /* -maskpolyatails */
  callinfo->simfilterparam.maskpolyAtails = false;
  if (!gthconsensus_parsing) {
    optmaskpolyatails = gt_option_new_bool("maskpolyatails", "mask poly(A) "
                                           "tails in cDNA/EST files",
                                           &callinfo->simfilterparam
                                           .maskpolyAtails,
                                           DEFAULT_MASKPOLYATAILS);
    gt_option_is_extended_option(optmaskpolyatails);
    gt_option_parser_add_option(op, optmaskpolyatails);
  }

  /* -proteinsmap */
  optproteinsmap = gt_option_new_string(PROTEINSMAP_OPT_CSTR, "specify smap "
                                        "file used for protein files",
                                        gth_input_proteinsmap(input),
                                        DEFAULT_PROTEINSMAP);
  gt_option_is_extended_option(optproteinsmap);
  gt_option_parser_add_option(op, optproteinsmap);

  /* -noautoindex */
  if (!gthconsensus_parsing) {
    optnoautoindex = gt_option_new_bool("noautoindex", "do not create indices "
                                     "automatically\nexcept for the .dna.* "
                                     "files used for the DP.\nexistence is not "
                                     "tested before an index is actually used!",
                                     &callinfo->simfilterparam.noautoindex,
                                     DEFAULT_NOAUTOINDEX);
    gt_option_is_extended_option(optnoautoindex);
    gt_option_parser_add_option(op, optnoautoindex);
  }

  /* -createindicesonly */
  optcreateindicesonly = gt_option_new_bool("createindicesonly", "stop program "
                                         "flow after the indices have been "
                                         "created",
                                         &callinfo->simfilterparam
                                         .createindicesonly,
                                         DEFAULT_CREATEINDICESONLY);
  gt_option_is_extended_option(optcreateindicesonly);
  gt_option_parser_add_option(op, optcreateindicesonly);

  /* -skipindexcheck */
  optskipindexcheck = gt_option_new_bool("skipindexcheck", "skip index check "
                                         "(in preprocessing phase)",
                                         &callinfo->simfilterparam
                                         .skipindexcheck,
                                         DEFAULT_SKIPINDEXCHECK);
  gt_option_is_extended_option(optskipindexcheck);
  gt_option_parser_add_option(op, optskipindexcheck);

  /* -minmatchlen */
  if (!gthconsensus_parsing) {
    optminmatchlen = gt_option_new_ulong_min(MINMATCHLEN_OPT_CSTR, "specify "
                                          "minimum match length (cDNA "
                                          "matching)", &callinfo->simfilterparam
                                          .minmatchlength,
                                          DEFAULT_MINMATCHLENGTH, 1);
    gt_option_parser_add_option(op, optminmatchlen);
  }

  /* -seedlength */
  if (!gthconsensus_parsing) {
    optseedlength = gt_option_new_ulong_min(SEEDLENGTH_OPT_CSTR, "specify the "
                                         "seed length (cDNA matching)",
                                         &callinfo->simfilterparam.seedlength,
                                         DEFAULT_SEEDLENGTH, 1);
    gt_option_parser_add_option(op, optseedlength);
  }

  /* -exdrop */
  if (!gthconsensus_parsing) {
    optexdrop = gt_option_new_ulong_min("exdrop", "specify the Xdrop value for "
                                     "edit distance extension (cDNA matching)",
                                     &callinfo->simfilterparam.exdrop,
                                     DEFAULT_EXDROP, 1);
    gt_option_parser_add_option(op, optexdrop);
  }

  /* -prminmatchlen */
  if (!gthconsensus_parsing) {
    optprminmatchlen = gt_option_new_ulong_min(PRMINMATCHLEN_OPT_CSTR,
                                               "specify minimum match length "
                                               "(protein matches)",
                                               &callinfo->simfilterparam
                                               .prminmatchlen,
                                               DEFAULT_PRMINMATCHLEN, 1);
    gt_option_parser_add_option(op, optprminmatchlen);
  }

  /* -prseedlength */
  if (!gthconsensus_parsing) {
    optprseedlength = gt_option_new_ulong_min(PRSEEDLENGTH_OPT_CSTR, "specify "
                                           "seed length (protein matching)",
                                           &callinfo->simfilterparam
                                           .prseedlength, DEFAULT_PRSEEDLENGTH,
                                           1);
    gt_option_parser_add_option(op, optprseedlength);
  }

  /* -prhdist */
  if (!gthconsensus_parsing) {
    optprhdist = gt_option_new_ulong_min("prhdist", "specify Hamming distance "
                                      "(protein matching)",
                                      &callinfo->simfilterparam.prhdist,
                                      DEFAULT_PRHDIST, 1);
    gt_option_parser_add_option(op, optprhdist);
  }

  /* -online */
  if (!gthconsensus_parsing) {
    optonline = gt_option_new_bool("online", "run the similarity filter online "
                                "without using the complete index (increases "
                                "runtime)", &callinfo->simfilterparam.online,
                                DEFAULT_ONLINE);
    gt_option_is_extended_option(optonline);
    gt_option_parser_add_option(op, optonline);
  }

  /* -inverse */
  if (!gthconsensus_parsing) {
    optinverse = gt_option_new_bool("inverse", "invert query and index in "
                                    "vmatch call",
                                    &callinfo->simfilterparam.inverse,
                                    DEFAULT_INVERSE);
    gt_option_is_extended_option(optinverse);
    gt_option_parser_add_option(op, optinverse);
  }

  /* -exact */
  if (!gthconsensus_parsing) {
    optexact = gt_option_new_bool("exact", "use exact matches in the "
                                  "similarity filter",
                                  &callinfo->simfilterparam.exact,
                                  DEFAULT_EXACT);
    gt_option_is_extended_option(optexact);
    gt_option_parser_add_option(op, optexact);
  }

  /* -edist */
  if (!gthconsensus_parsing) {
    optedist = gt_option_new_bool("edist", "use edist instead of exdrop in "
                               "simfilter", &callinfo->simfilterparam.edist,
                               DEFAULT_EDIST);
    gt_option_is_extended_option(optedist);
    gt_option_is_development_option(optedist);
    gt_option_parser_add_option(op, optedist);
  }

  /* -maxnumofmatches */
  if (!gthconsensus_parsing) {
    optmaxnumofmatches = gt_option_new_ulong("maxnumofmatches", "set the "
                                             "maximum number of matches per "
                                             "genomic file. I.e., only that "
                                             "number of matches will be saved "
                                             "for every reference sequence",
                                             &callinfo->simfilterparam
                                             .maxnumofmatches,
                                             DEFAULT_MAXNUMOFMATCHES);
    gt_option_is_extended_option(optmaxnumofmatches);
    gt_option_is_development_option(optmaxnumofmatches);
    gt_option_parser_add_option(op, optmaxnumofmatches);
  }

  /* -fragweightfactor */
  if (!gthconsensus_parsing) {
    optfragweightfactor = gt_option_new_double_min("fragweightfactor", "set "
                                                   "the weight factor with "
                                                   "which the fragments are "
                                                   "multiplied before the "
                                                   "global chaining is "
                                                   "performed",
                                                   &callinfo->fragweightfactor,
                                                   DEFAULT_FRAGWEIGHTFACTOR,
                                                   0.1);
    gt_option_is_development_option(optfragweightfactor);
    gt_option_parser_add_option(op, optfragweightfactor);
  }

  /* -gcmaxgapwidth */
  if (!gthconsensus_parsing) {
    optgcmaxgapwidth = gt_option_new_uint("gcmaxgapwidth", "set the maximum "
                                          "gap width for global chains\n"
                                          "defines approximately the maximum "
                                          "intron length\nset to 0 to allow "
                                          "for unlimited length\nin order to "
                                          "avoid false-positive exons (lonely "
                                          "exons) at the sequence ends, it is "
                                          "very important to set this "
                                          "parameter appropriately!",
                                          &callinfo->gcmaxgapwidth,
                                          DEFAULT_GCMAXGAPWIDTH);
    gt_option_parser_add_option(op, optgcmaxgapwidth);
  }

  /* -rare */
  if (!gthconsensus_parsing) {
    optrare = gt_option_new_ulong("rare", "set the maximum number of (rare) "
                               "matches. I.e., for any given start position in "
                               "the reference sequence, maximally ``rare'' "
                               "number of matches are allowed.",
                               &callinfo->simfilterparam.rare, DEFAULT_RARE);
    gt_option_is_development_option(optrare);
    gt_option_parser_add_option(op, optrare);
  }

  /* -gcmincoverage */
  if (!gthconsensus_parsing) {
    optgcmincoverage = gt_option_new_uint_max("gcmincoverage", "set the "
                                              "minimum coverage of global "
                                              "chains regarding to the "
                                              "reference sequence",
                                              &callinfo->gcmincoverage,
                                              DEFAULT_GCMINCOVERAGE, 100);
    gt_option_parser_add_option(op, optgcmincoverage);
  }

  /* -paralogs */
  if (!gthconsensus_parsing) {
    optparalogs = gt_option_new_bool("paralogs", "compute paralogous genes "
                                  "(different chaining procedure)",
                                  &callinfo->simfilterparam.paralogs,
                                  DEFAULT_PARALOGS);
    gt_option_parser_add_option(op, optparalogs);
  }

  /* -enrichchains */
  if (!gthconsensus_parsing) {
    optenrichchains = gt_option_new_bool("enrichchains", "enrich genomic "
                                         "sequence part of global chains with "
                                         "additional matches",
                                         &callinfo->simfilterparam.enrichchains,
                                         DEFAULT_ENRICHCHAINS);
    gt_option_is_development_option(optenrichchains);
    gt_option_parser_add_option(op, optenrichchains);
  }

  /* -stopafterchaining */
  if (!gthconsensus_parsing) {
    optstopafterchaining = gt_option_new_bool("stopafterchaining", "stop gth "
                                           "after chaining phase and show "
                                           "stored global chains",
                                           &callinfo->simfilterparam
                                           .stopafterchaining,
                                           DEFAULT_STOPAFTERCHAINING);
    gt_option_is_development_option(optstopafterchaining);
    gt_option_parser_add_option(op, optstopafterchaining);
  }

  /* -introncutout */
  if (!gthconsensus_parsing) {
    optintroncutout = gt_option_new_bool("introncutout", "enable the intron "
                                      "cutout technique",
                                      &callinfo->simfilterparam.introncutoutinfo
                                      .introncutout, DEFAULT_INTRONCUTOUT);
    gt_option_parser_add_option(op, optintroncutout);
  }

  /* -autointroncutout */
  if (!gthconsensus_parsing) {
    optautointroncutout = gt_option_new_uint("autointroncutout", "set the "
                                          "automatic intron cutout matrix size "
                                          "in megabytes and enable the "
                                          "automatic intron cutout technique",
                                           &callinfo->simfilterparam
                                           .introncutoutinfo
                                           .autoicmaxmatrixsize,
                                           DEFAULT_AUTOICMAXMATRIXSIZE);
    gt_option_parser_add_option(op, optautointroncutout);
  }

  /* -icinitialdelta */
  if (!gthconsensus_parsing) {
    opticinitialdelta = gt_option_new_uint(ICINITIALDELTA_OPT_CSTR, "set the "
                                        "initial delta used for intron cutouts",
                                        &callinfo->simfilterparam
                                        .introncutoutinfo.icinitialdelta,
                                        DEFAULT_ICINITIALDELTA);
    gt_option_is_extended_option(opticinitialdelta);
    gt_option_parser_add_option(op, opticinitialdelta);
  }

  /* -iciterations */
  if (!gthconsensus_parsing) {
    opticiterations = gt_option_new_uint("iciterations", "set the number of "
                                      "intron cutout iterations",
                                      &callinfo->simfilterparam.introncutoutinfo
                                      .iciterations, DEFAULT_ICITERATIONS);
    gt_option_is_extended_option(opticiterations);
    gt_option_parser_add_option(op, opticiterations);
  }

  /* -icdeltaincrease */
  if (!gthconsensus_parsing) {
    opticdeltaincrease = gt_option_new_uint("icdeltaincrease", "set the delta "
                                         "increase during every iteration",
                                         &callinfo->simfilterparam
                                         .introncutoutinfo.icdeltaincrease,
                                         DEFAULT_ICDELTAINCREASE);
    gt_option_is_extended_option(opticdeltaincrease);
    gt_option_parser_add_option(op, opticdeltaincrease);
  }

  /* -icminremintronlen */
  if (!gthconsensus_parsing) {
    opticminremlength = gt_option_new_uint("icminremintronlen", "set the "
                                           "minimum remaining intron length "
                                           "for an intron to be cut out",
                                           &callinfo->simfilterparam
                                           .introncutoutinfo
                                           .icminremintronlength,
                                           DEFAULT_ICMINREMLENGTH);
    gt_option_is_extended_option(opticminremlength);
    gt_option_parser_add_option(op, opticminremlength);
  }

  /* -nou12intronmodel */
  if (!gthconsensus_parsing) {
    optnou12intronmodel = gt_option_new_bool("nou12intronmodel", "disable the "
                                          "U12-type intron model",
                                          &nou12intronmodel,
                                          DEFAULT_DISABLEU12INTRONMODEL);
    gt_option_is_extended_option(optnou12intronmodel);
    gt_option_parser_add_option(op, optnou12intronmodel);
  }

  /* -u12donorprob */
  if (!gthconsensus_parsing) {
    optu12donorprob = gt_option_new_probability("u12donorprob", "set the "
                                                "probability for perfect "
                                                "U12-type donor sites",
                                                &u12donorprob,
                                                DEFAULT_U12_TYPEDONORPROB);
    gt_option_is_extended_option(optu12donorprob);
    gt_option_parser_add_option(op, optu12donorprob);
  }

  /* -u12donorprob1mism */
  if (!gthconsensus_parsing) {
    optu12donorprob1mism = gt_option_new_probability("u12donorprob1mism", "set "
                                                  "the prob. for U12-type "
                                                  "donor w. 1 mismatch",
                                                  &u12donorprob1mism,
                                          DEFAULT_U12_TYPEDONORPROBONEMISMATCH);
    gt_option_is_extended_option(optu12donorprob1mism);
    gt_option_parser_add_option(op, optu12donorprob1mism);
  }

  /* -probies */
  if (!gthconsensus_parsing) {
    optprobies = gt_option_new_probability("probies", "set the initial exon "
                                           "state probability",
                                           &callinfo->dp_options_est->probies,
                                           GTH_DEFAULT_PROBIES);
    gt_option_is_extended_option(optprobies);
    gt_option_parser_add_option(op, optprobies);
  }

  /* -probdelgen */
  if (!gthconsensus_parsing) {
    optprobdelgen = gt_option_new_probability("probdelgen", "set the genomic "
                                              "sequence deletion probability",
                                              &callinfo->dp_options_est
                                              ->probdelgen,
                                              GTH_DEFAULT_PROBDELGEN);
    gt_option_is_extended_option(optprobdelgen);
    gt_option_parser_add_option(op, optprobdelgen);
  }

  /* -identityweight */
  if (!gthconsensus_parsing) {
    optidentityweight = gt_option_new_double("identityweight", "set the pairs "
                                             "of identical characters weight",
                                             &callinfo->dp_options_est
                                             ->identityweight,
                                             GTH_DEFAULT_IDENTITYWEIGHT);
    gt_option_is_extended_option(optidentityweight);
    gt_option_parser_add_option(op, optidentityweight);
  }

  /* -mismatchweight */
  if (!gthconsensus_parsing) {
    optmismatchweight = gt_option_new_double("mismatchweight", "set the weight "
                                             "for mismatching characters",
                                             &callinfo->dp_options_est
                                             ->mismatchweight,
                                             GTH_DEFAULT_MISMATCHWEIGHT);
    gt_option_is_extended_option(optmismatchweight);
    gt_option_parser_add_option(op, optmismatchweight);
  }

  /* -undetcharweight */
  if (!gthconsensus_parsing) {
    optundetcharweight = gt_option_new_double("undetcharweight", "set the "
                                              "weight for undetermined "
                                              "characters",
                                              &callinfo->dp_options_est
                                              ->undetcharweight,
                                              GTH_DEFAULT_UNDETCHARWEIGHT);
    gt_option_is_extended_option(optundetcharweight);
    gt_option_parser_add_option(op, optundetcharweight);
  }

  /* -deletionweight */
  if (!gthconsensus_parsing) {
    optdeletionweight = gt_option_new_double("deletionweight", "set the weight "
                                             "for deletions",
                                             &callinfo->dp_options_est
                                             ->deletionweight,
                                             GTH_DEFAULT_DELETIONWEIGHT);
    gt_option_is_extended_option(optdeletionweight);
    gt_option_parser_add_option(op, optdeletionweight);
  }

  /* -freeintrontrans */
  optfreeintrontrans = gt_option_new_bool("freeintrontrans", "free transitions "
                                          "between intron states",
                                          &callinfo->dp_options_core
                                          ->freeintrontrans,
                                          GTH_DEFAULT_FREEINTRONTRANS);
  gt_option_is_development_option(optfreeintrontrans);
  gt_option_parser_add_option(op, optfreeintrontrans);

  /* -dpminexonlen */
  if (!gthconsensus_parsing) {
    optdpminexonlength = gt_option_new_uint_min("dpminexonlen", "set the "
                                                "minimum exon length for the "
                                                "DP", &callinfo->dp_options_core
                                                      ->dpminexonlength,
                                                GTH_DEFAULT_DPMINEXONLENGTH, 1);
    gt_option_is_extended_option(optdpminexonlength);
    gt_option_parser_add_option(op, optdpminexonlength);
  }

  /* -dpminintronlen */
  if (!gthconsensus_parsing) {
    optdpminintronlength = gt_option_new_uint_min(DPMININTRONLENGTH_OPT_CSTR,
                                               "set the minimum intron length "
                                               "for the DP",
                                               &callinfo->dp_options_core
                                               ->dpminintronlength,
                                               GTH_DEFAULT_DPMININTRONLENGTH,
                                               1);
    gt_option_is_extended_option(optdpminintronlength);
    gt_option_parser_add_option(op, optdpminintronlength);
  }

  /* -shortexonpenal */
  if (!gthconsensus_parsing) {
    optshortexonpenalty = gt_option_new_double_min("shortexonpenal", "set the "
                                                "short exon penalty",
                                                &callinfo->dp_options_core
                                                ->shortexonpenalty,
                                                GTH_DEFAULT_SHORTEXONPENALTY,
                                                0.0);
    gt_option_is_extended_option(optshortexonpenalty);
    gt_option_parser_add_option(op, optshortexonpenalty);
  }

  /* -shortintronpenal */
  if (!gthconsensus_parsing) {
    optshortintronpenalty =
      gt_option_new_double_min("shortintronpenal",
                               "set the short intron penalty",
                               &callinfo->dp_options_core->shortintronpenalty,
                               GTH_DEFAULT_SHORTINTRONPENALTY, 0.0);
    gt_option_is_extended_option(optshortintronpenalty);
    gt_option_parser_add_option(op, optshortintronpenalty);
  }

  /* -wzerotransition */
  if (!gthconsensus_parsing) {
    optwzerotransition = gt_option_new_uint("wzerotransition", "set the zero "
                                         "transition weights window size",
                                         &callinfo->dp_options_est
                                         ->wzerotransition,
                                         GTH_DEFAULT_WZEROTRANSITION);
    gt_option_is_extended_option(optwzerotransition);
    gt_option_parser_add_option(op, optwzerotransition);
  }

  /* -wdecreasedoutput */
  if (!gthconsensus_parsing) {
    optwdecreasedoutput = gt_option_new_uint("wdecreasedoutput", "set the "
                                          "decreased output weights window "
                                          "size", &callinfo->dp_options_est
                                          ->wdecreasedoutput,
                                          GTH_DEFAULT_WDECREASEDOUTPUT);
    gt_option_is_extended_option(optwdecreasedoutput);
    gt_option_parser_add_option(op, optwdecreasedoutput);
  }

  /* -detectsmallexons */
  if (!gthconsensus_parsing) {
    optdetectsmallexons = gt_option_new_bool("detectsmallexons", "detect small "
                                            "exons with additional dynamic "
                                            "programming calls",
                                            &callinfo->dp_options_est
                                            ->detectsmallexons,
                                            GTH_DEFAULT_DETECTSMALLEXONS);
    gt_option_is_development_option(optdetectsmallexons);
    gt_option_parser_add_option(op, optdetectsmallexons);
  }

  /* -noicinintroncheck */
  if (!gthconsensus_parsing) {
    optnoicinintroncheck = gt_option_new_bool("noicinintroncheck", "disable "
                                           "intron cutout in intron check\nin "
                                           "such cases include exon consisting "
                                           "of deletions instead",
                                           &callinfo->dp_options_core
                                           ->noicinintroncheck,
                                           GTH_DEFAULT_NOICININTRONCHECK);
    gt_option_is_extended_option(optnoicinintroncheck);
    gt_option_is_development_option(optnoicinintroncheck);
    gt_option_parser_add_option(op, optnoicinintroncheck);
  }

  /* -proteinexonpenal */
  if (!gthconsensus_parsing) {
    optproteinexonpenal = gt_option_new_bool("proteinexonpenal", "use short "
                                             "exon penalties in protein DP",
                                             &callinfo->proteinexonpenal,
                                             DEFAULT_PROTEINEXONPENAL);
    gt_option_is_extended_option(optproteinexonpenal);
    gt_option_is_development_option(optproteinexonpenal);
    gt_option_parser_add_option(op, optproteinexonpenal);
  }

  /* -leadcutoffsmode */
  if (!gthconsensus_parsing) {
    optleadcutoffsmode = gt_option_new_choice("leadcutoffsmode", "set the "
                                              "cutoffs mode for leading "
                                              "bases\ncan be either RELAXED, "
                                              "STRICT, or MINIMAL",
                                              leadcutoffsmode, cutoff_modes[0],
                                              cutoff_modes);
    gt_option_is_extended_option(optleadcutoffsmode);
    gt_option_parser_add_option(op, optleadcutoffsmode);
  }

  /* -termcutoffsmode */
  if (!gthconsensus_parsing) {
    opttermcutoffsmode = gt_option_new_choice("termcutoffsmode", "set the "
                                              "cutoffs mode for terminal "
                                              "bases\ncan be either RELAXED, "
                                              "STRICT, or MINIMAL",
                                              termcutoffsmode, cutoff_modes[1],
                                              cutoff_modes);
    gt_option_is_extended_option(opttermcutoffsmode);
    gt_option_parser_add_option(op, opttermcutoffsmode);
  }

  /* -cutoffsminexonlen */
  if (!gthconsensus_parsing) {
    optcutoffsminexonlength = gt_option_new_uint_min("cutoffsminexonlen",
                                                  "set the cutoffs minimum "
                                                  "exon length",
                                                  &callinfo->dp_options_postpro
                                                  ->cutoffsminexonlen,
                                                  GTH_DEFAULT_CUTOFFSMINEXONLEN,
                                                  1);
    gt_option_is_extended_option(optcutoffsminexonlength);
    gt_option_parser_add_option(op, optcutoffsminexonlength);
  }

  /* -scoreminexonlen */
  if (!gthconsensus_parsing) {
    optscoreminexonlength = gt_option_new_uint_min("scoreminexonlen", "set the "
                                                "score minimum exon length",
                                                &callinfo->dp_options_postpro
                                                ->scoreminexonlen,
                                                GTH_DEFAULT_SCOREMINEXONLEN, 1);
    gt_option_is_extended_option(optscoreminexonlength);
    gt_option_parser_add_option(op, optscoreminexonlength);
  }

  /* -minaveragessp */
  if (!gthconsensus_parsing) {
    optminaveragessp = gt_option_new_probability("minaveragessp", "set the "
                                              "minimum average splice site "
                                              "prob.", &callinfo->minaveragessp,
                                              DEFAULT_MINAVERAGESSP);
    gt_option_is_extended_option(optminaveragessp);
    gt_option_parser_add_option(op, optminaveragessp);
  }

  /* add spliced alignment filter options */
  gth_sa_filter_register_options(op, callinfo->sa_filter,
                                 !gthconsensus_parsing);

  /* -intermediate */
  optintermediate = gt_option_new_bool("intermediate", "stop after calculation "
                                       "of spliced alignments and output "
                                       "results in reusable XML format. Do not "
                                       "process this output yourself, use the "
                                       "``normal'' XML output instead!",
                                       &callinfo->intermediate,
                                       DEFAULT_INTERMEDIATE);
  gt_option_parser_add_option(op, optintermediate);

  /* -sortags */
  optsortags = gt_option_new_bool("sortags", "sort alternative gene structures "
                               "according to the weighted mean of the average "
                               "exon score and the average splice site "
                               "probability", &callinfo->out->sortags,
                               DEFAULT_SORTAGS);
  gt_option_is_extended_option(optsortags);
  gt_option_parser_add_option(op, optsortags);

  /* -sortagswf */
  optsortagswf = gt_option_new_double_min("sortagswf", "set the weight factor "
                                          "for the sorting of AGSs",
                                          &callinfo->out->sortagswf,
                                          DEFAULT_SORTAGSWF, 0.001);
  gt_option_is_extended_option(optsortagswf);
  gt_option_parser_add_option(op, optsortagswf);

  /* -disableclustersas */
  optdisableclustersas = gt_option_new_bool("disableclustersas", "disable the "
                                         "clustering of spliced alignments in "
                                         "the consensus phase (should not "
                                         "change the results)",
                                         &callinfo->disableclustersas,
                                         DEFAULT_DISABLECLUSTERSAS);
  gt_option_is_extended_option(optdisableclustersas);
  gt_option_is_development_option(optdisableclustersas);
  gt_option_parser_add_option(op, optdisableclustersas);

  /* -cdnaforwardonly */
  if (!gthconsensus_parsing) {
    optcdnaforwardonly = gt_option_new_bool("cdnaforwardonly", "consider only "
                                         "forward strand of cDNAs",
                                         &callinfo->cdnaforwardonly, false);
    gt_option_is_development_option(optcdnaforwardonly);
    gt_option_parser_add_option(op, optcdnaforwardonly);
  }

  /* -exondistri */
  optexondistri = gt_option_new_bool("exondistri", "show the exon length "
                                  "distribution", &exondistri, false);
  gt_option_is_extended_option(optexondistri);
  gt_option_parser_add_option(op, optexondistri);

  /* -introndistri */
  optintrondistri = gt_option_new_bool("introndistri", "show the intron length "
                                    "distribution", &introndistri, false);
  gt_option_is_extended_option(optintrondistri);
  gt_option_parser_add_option(op, optintrondistri);

  /* -refseqcovdistri */
  if (!gthconsensus_parsing) {
    optrefseqcovdistri = gt_option_new_bool("refseqcovdistri", "show the "
                                            "reference sequence coverage "
                                            "distribution", &refseqcovdistri,
                                            false);
    gt_option_is_extended_option(optrefseqcovdistri);
    gt_option_parser_add_option(op, optrefseqcovdistri);
  }

  /* -matchnumdistri */
  if (!gthconsensus_parsing) {
    optmatchnumdistri = gt_option_new_bool("matchnumdistri", "show the "
                                        "distribution of matches\n(per genomic "
                                        "file and per reference sequence)",
                                        &matchnumdistri, false);
    gt_option_is_extended_option(optmatchnumdistri);
    gt_option_is_development_option(optmatchnumdistri);
    gt_option_parser_add_option(op, optmatchnumdistri);
  }

  /* -first */
  if (!gthconsensus_parsing) {
    optfirstalshown = gt_option_new_uint("first", "set the maximum number of "
                                      "spliced alignments per genomic DNA "
                                      "input. Set to 0 for unlimited number.",
                                      &callinfo->firstalshown,
                                      DEFAULT_FIRSTALSHOWN);
    gt_option_parser_add_option(op, optfirstalshown);
  }

  /* -showeops */
  if (!gthconsensus_parsing) {
    optshoweops = gt_option_new_bool("showeops", "show complete array of multi "
                                  "edit operations after each DP",
                                  &callinfo->out->showeops, DEFAULT_SHOWEOPS);
    gt_option_is_extended_option(optshoweops);
    gt_option_is_development_option(optshoweops);
    gt_option_parser_add_option(op, optshoweops);
  }

  /* mandatory options  */
  if (!gthconsensus_parsing) {
    gt_option_is_mandatory(optgenomic);
    gt_option_is_mandatory_either(optcdna, optprotein);
  }

  /* option exclusions */
  if (optforward && optreverse)
    gt_option_exclude(optforward, optreverse);
  if (optspecies && optbssm)
    gt_option_exclude(optspecies, optbssm);
  if (optmincutoffs && optleadcutoffsmode)
    gt_option_exclude(optmincutoffs, optleadcutoffsmode);
  if (optmincutoffs && opttermcutoffsmode)
    gt_option_exclude(optmincutoffs, opttermcutoffsmode);
  if (optexact && optseedlength)
    gt_option_exclude(optexact, optseedlength);
  if (optexact && optexdrop)
    gt_option_exclude(optexact, optexdrop);
  if (optgs2out && optshowseqnums)
    gt_option_exclude(optgs2out, optshowseqnums);
  if (optxmlout && optshowseqnums)
    gt_option_exclude(optxmlout, optshowseqnums);
  if (optexact && optedist)
    gt_option_exclude(optexact, optedist);
  if (optintroncutout && optautointroncutout)
    gt_option_exclude(optintroncutout, optautointroncutout);
  if (optmaskpolyatails && optnoautoindex)
    gt_option_exclude(optmaskpolyatails, optnoautoindex);
  if (optproteinsmap && optnoautoindex)
    gt_option_exclude(optproteinsmap, optnoautoindex);
  if (optfrompos && optonline)
    gt_option_exclude(optfrompos, optonline);
  if (optnou12intronmodel && optu12donorprob)
    gt_option_exclude(optnou12intronmodel, optu12donorprob);
  if (optnou12intronmodel && optu12donorprob1mism)
    gt_option_exclude(optnou12intronmodel, optu12donorprob1mism);
  if (optnoautoindex && optcreateindicesonly)
    gt_option_exclude(optnoautoindex, optcreateindicesonly);
  if (optnoautoindex && optskipindexcheck)
    gt_option_exclude(optnoautoindex, optskipindexcheck);
  if (optcreateindicesonly && optskipindexcheck)
    gt_option_exclude(optcreateindicesonly, optskipindexcheck);
  if (optskipalignmentout && optintermediate)
    gt_option_exclude(optskipalignmentout, optintermediate);
  if (optxmlout && optgff3out)
    gt_option_exclude(optxmlout, optgff3out);

  /* option implications (single) */
  if (opttopos && optfrompos)
    gt_option_imply(opttopos, optfrompos);
  if (optwidth && optfrompos)
    gt_option_imply(optwidth, optfrompos);
  if (optverboseseqs && optverbose)
    gt_option_imply(optverboseseqs, optverbose);
  if (optsortagswf && optsortags)
    gt_option_imply(optsortagswf, optsortags);

  /* option implications (either 2) */
  if (optgff3out && optskipalignmentout && optintermediate)
    gt_option_imply_either_2(optgff3out, optskipalignmentout, optintermediate);
  if (opticinitialdelta && optintroncutout && optautointroncutout) {
    gt_option_imply_either_2(opticinitialdelta, optintroncutout,
                          optautointroncutout);
  }
  if (opticiterations && optintroncutout && optautointroncutout) {
    gt_option_imply_either_2(opticiterations, optintroncutout,
                          optautointroncutout);
  }
  if (opticdeltaincrease && optintroncutout && optautointroncutout) {
    gt_option_imply_either_2(opticdeltaincrease, optintroncutout,
                          optautointroncutout);
  }
  if (opticminremlength && optintroncutout && optautointroncutout) {
    gt_option_imply_either_2(opticminremlength, optintroncutout,
                          optautointroncutout);
  }
  if (optintermediate && optxmlout && optgff3out) {
    gt_option_imply_either_2(optintermediate, optxmlout, optgff3out);
  }

  /* set mail addresse */
  gt_option_parser_set_mailaddress(op, "<gremme@zbh.uni-hamburg.de>");

  /* set comment function */
  gt_option_parser_set_comment_func(op, show_gth_help_trailer, NULL);

  /* parse options */
  if (gthconsensus_parsing) {
    oprval = gt_option_parser_parse(op, parsed_args, argc, argv, show_version,
                                    err);
  }
  else {
    gt_option_parser_set_max_args(op, 0);
    oprval = gt_option_parser_parse(op, parsed_args, argc, argv, show_version,
                                    err);
  }

  if (oprval == OPTIONPARSER_OK &&
      callinfo->out->showintronmaxlen > 0 &&
      callinfo->out->showintronmaxlen < DOUBLEALIGNMENTLINEWIDTH) {
     gt_error_set(err, "argument to option -%s must be 0 or integer larger "
                       "than %u", SHOWINTRONMAXLEN_OPT_CSTR,
                  DOUBLEALIGNMENTLINEWIDTH);
     oprval = OPTIONPARSER_ERROR;
  }

  if (oprval == OPTIONPARSER_OK) {
    /* process genomic files */
    for (i = 0; i < gt_str_array_size(genomic_files); i++)
      gth_input_add_genomic_file(input, gt_str_array_get(genomic_files, i));

    /* process cdna files */
    for (i = 0; i < gt_str_array_size(cdna_files); i++)
      gth_input_add_cdna_file(input, gt_str_array_get(cdna_files, i));

    /* process protein files */
    for (i = 0; i < gt_str_array_size(protein_files); i++)
      gth_input_add_protein_file(input, gt_str_array_get(protein_files, i));
  }

  if (oprval == OPTIONPARSER_OK && gt_str_length(parsed_species)) {
    ret = findspeciesnum(gt_str_get(parsed_species), err);
    if (ret == -1)
      oprval = OPTIONPARSER_ERROR;
    if (oprval == OPTIONPARSER_OK) {
      callinfo->speciesnum = ret;
      /* also setting the bssm parameter file */
      gt_assert(!gt_str_length(gth_input_bssmfile(input)));
      gt_str_append_cstr(gth_input_bssmfile(input), speciestab[ret]);
    }
  }

  /* check translation table number */
  /* XXX: the maximal translation tables is set to 15, because that is the
     current limit in Vmatch, this should be changed in Vmatch and the check
     refactored into the GtTranslator class */
  if (oprval == OPTIONPARSER_OK &&
      (callinfo->translationtable < 1  || callinfo->translationtable > 15 ||
       callinfo->translationtable == 7 || callinfo->translationtable == 8)) {
    gt_error_set(err, "wrong argument '%u' to option \"-%s\": "
                 "must be number in the range [1,15] except for 7 and 8",
                 callinfo->translationtable,
                 gt_option_get_name(opttranslationtable));
    oprval = OPTIONPARSER_ERROR;
  }

  if (!gthconsensus_parsing) {
    /* some checking necessary if -frompos is set */
    if (oprval == OPTIONPARSER_OK && gt_option_is_set(optfrompos)) {
      /* if -topos and -length are used at the same time throw an error */
      if (gt_option_is_set(opttopos) && gt_option_is_set(optwidth)) {
        gt_error_set(err, "-%s and -%s can not be used simultaneously",
                  TOPOS_OPT_CSTR, WIDTH_OPT_CSTR);
        oprval = OPTIONPARSER_ERROR;
      }
      /* if neither -topos nor -length are used throw an error */
      if (oprval == OPTIONPARSER_OK && !gt_option_is_set(opttopos) &&
          !gt_option_is_set(optwidth)) {
        gt_error_set(err, "if -%s is used either -%s or -%s have to be used",
                  FROMPOS_OPT_CSTR, TOPOS_OPT_CSTR, WIDTH_OPT_CSTR);
        oprval = OPTIONPARSER_ERROR;
      }
    }

    /* make sure we have only one genomic file */
    if (oprval == OPTIONPARSER_OK && gt_option_is_set(optfrompos) &&
        gth_input_num_of_gen_files(input) > 1) {
      gt_error_set(err, "-%s only allowed if exactly 1 genomic file is "
                        "specified", FROMPOS_OPT_CSTR);
      oprval = OPTIONPARSER_ERROR;
    }

    /* if optfrompos is set, we set optinverse automatically. This is more
       intuitive for the user (instead of requiring him to set -inverse) */
    if (oprval == OPTIONPARSER_OK && gt_option_is_set(optfrompos) &&
        !gt_option_is_set(optinverse)) {
      callinfo->simfilterparam.inverse = true;
    }
  }

  /* check necessary for intron cutout technique */
  if (oprval == OPTIONPARSER_OK && !gthconsensus_parsing &&
      2 * callinfo->simfilterparam.introncutoutinfo.icinitialdelta <
      callinfo->dp_options_core->dpminintronlength) {
    gt_error_set(err, "2 * %s(=%u) must be >= %s(=%u)",
              ICINITIALDELTA_OPT_CSTR,
              callinfo->simfilterparam.introncutoutinfo.icinitialdelta,
              DPMININTRONLENGTH_OPT_CSTR,
              callinfo->dp_options_core->dpminintronlength);
    oprval = OPTIONPARSER_ERROR;
  }

  /* check number of genomic files */
  if (oprval == OPTIONPARSER_OK &&
      gth_input_num_of_gen_files(input) >= VM_MAXNUMBEROFFILES) {
    gt_error_set(err, "maximal number of genomic files is %u",
                 VM_MAXNUMBEROFFILES);
    oprval = OPTIONPARSER_ERROR;
  }

  /* make sure minmatchlen is >= seedlength */
  if (oprval == OPTIONPARSER_OK && !gthconsensus_parsing &&
      callinfo->simfilterparam.minmatchlength <
      callinfo->simfilterparam.seedlength) {
    gt_error_set(err, "-%s %lu must be >= -%s %lu", MINMATCHLEN_OPT_CSTR,
              callinfo->simfilterparam.minmatchlength, SEEDLENGTH_OPT_CSTR,
              callinfo->simfilterparam.seedlength);
    oprval = OPTIONPARSER_ERROR;
  }

  /* make sure prminmatchlen is >= prseedlength */
  if (oprval == OPTIONPARSER_OK && !gthconsensus_parsing &&
      callinfo->simfilterparam.prminmatchlen <
      callinfo->simfilterparam.prseedlength) {
    gt_error_set(err, "-%s %lu must be >= -%s %lu", PRMINMATCHLEN_OPT_CSTR,
              callinfo->simfilterparam.prminmatchlen, PRSEEDLENGTH_OPT_CSTR,
              callinfo->simfilterparam.prseedlength);
    oprval = OPTIONPARSER_ERROR;
  }

  /* check that given proteinsmap is not too long (is used as suffix for
     protein index files) */
  if (oprval == OPTIONPARSER_OK && optproteinsmap &&
      !gt_option_is_set(optproteinsmap)) {
    if (gt_str_length(gth_input_proteinsmap(input)) > MAXSUFFIXLEN) {
      gt_error_set(err, "\%s argument \"%s\" is too long (>MAXSUFFIXLEN=%u)",
                PROTEINSMAP_OPT_CSTR,
                gt_str_get(gth_input_proteinsmap(input)), MAXSUFFIXLEN);
      oprval = OPTIONPARSER_ERROR;
    }
  }

  if (oprval == OPTIONPARSER_OK) {
    gt_assert(!(forward && reverse));
    if (forward)
      gth_input_set_forward_only(input);
    if (reverse)
      gth_input_set_reverse_only(input);

    if (verbose) {
      callinfo->out->showverbose   = showverbose;
      callinfo->out->showverboseVM = showverboseVM;
    }

    if (!gthconsensus_parsing) {
      gth_splice_site_model_U12intronmodel_set_usage(callinfo
                                                     ->splice_site_model,
                                                     !nou12intronmodel);
      gth_splice_site_model_set_U12typedonorprob(callinfo->splice_site_model,
                                                 u12donorprob);
      gth_splice_site_model_set_U12typedonorprob_one_mismatch(callinfo
                                                            ->splice_site_model,
                                                             u12donorprob1mism);
    }

    if (exondistri)
      gth_stat_enable_exondistri(stat);
    if (introndistri)
      gth_stat_enable_introndistri(stat);
    if (refseqcovdistri)
      gth_stat_enable_refseqcovdistri(stat);
    if (matchnumdistri)
      gth_stat_enable_matchnumdistri(stat);
  }

  /* post-process cutoff mode options */
  if (!gthconsensus_parsing) {
    if (oprval == OPTIONPARSER_OK) {
      gt_assert(gt_str_length(leadcutoffsmode));
      mode = get_cutoffs_mode_from_table(gt_str_get(leadcutoffsmode));
      callinfo->dp_options_postpro->leadcutoffsmode = mode;
    }
    if (oprval == OPTIONPARSER_OK) {
      gt_assert(gt_str_length(termcutoffsmode));
      mode = get_cutoffs_mode_from_table(gt_str_get(termcutoffsmode));
      callinfo->dp_options_postpro->termcutoffsmode = mode;
    }
    if (oprval == OPTIONPARSER_OK && mincutoffs) {
      callinfo->dp_options_postpro->leadcutoffsmode = MINIMAL;
      callinfo->dp_options_postpro->termcutoffsmode = MINIMAL;
    }
  }

  /* assertions */
#ifndef NDEBUG
  if (oprval == OPTIONPARSER_OK && !gthconsensus_parsing) {
    /* at least one genomic file defined */
    gt_assert(gth_input_num_of_gen_files(input));
    /* at least one reference file defined */
    gt_assert(gth_input_num_of_ref_files(input));
  }
#endif

  /* save consensus files */
  if (oprval == OPTIONPARSER_OK && gthconsensus_parsing) {
    while (*parsed_args < argc) {
      gt_str_array_add_cstr(consensusfiles, argv[*parsed_args]);
      (*parsed_args)++;
    }
  }

  /* load bssm parameter */
  if (oprval == OPTIONPARSER_OK) {
    if (gt_str_length(gth_input_bssmfile(input))) {
      if (gth_splice_site_model_load_bssm(callinfo->splice_site_model,
                                          gt_str_get(gth_input_bssmfile(input)),
                                          err)) {
        oprval = OPTIONPARSER_ERROR;
      }
    }
  }

  if (OPTIONPARSER_OK) {
    gt_assert(*parsed_args == argc);
  }

  /* free */
  gt_str_delete(parsed_species);
  gt_str_delete(specbuf);
  gt_str_delete(termcutoffsmode);
  gt_str_delete(leadcutoffsmode);
  gt_str_delete(chaining_param);
  gt_str_array_delete(protein_files);
  gt_str_array_delete(cdna_files);
  gt_str_array_delete(genomic_files);
  gt_outputfileinfo_delete(ofi);
  gt_option_parser_delete(op);

  return oprval;
}
