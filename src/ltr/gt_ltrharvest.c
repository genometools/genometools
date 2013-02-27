/*
  Copyright (c) 2010-2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007      David Ellinghaus <d.ellinghaus@ikmb.uni-kiel.de>
  Copyright (c) 2007-2012 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "core/ma.h"
#include "core/option_api.h"
#include "core/unused_api.h"
#include "core/output_file_api.h"
#include "core/undef_api.h"
#include "core/versionfunc.h"
#include "extended/genome_node.h"
#include "extended/gff3_out_stream_api.h"
#include "ltr/ltr_four_char_motif.h"
#include "ltr/ltrharvest_stream.h"
#include "ltr/ltrharvest_fasta_out_stream.h"
#include "ltr/ltrharvest_tabout_stream.h"
#include "ltr/ltrharvest_tabout_visitor.h"
#include "ltr/gt_ltrharvest.h"

typedef struct {
  GtOutputFileInfo *ofi;
  GtFile *outfp;
  GtStr *str_indexname,
        *str_motif,
        *str_overlaps,
        *str_fastaoutputfilename,
        *str_fastaoutputfilenameinnerregion,
        *str_gff3filename;
  GtRange searchrange;
  unsigned long minseedlength,
                minltrlength,
                maxltrlength,
                mindistance,
                maxdistance,
                numofboundaries,
                offset;
  double similaritythreshold;
  int xdropbelowscore;
  GtXdropArbitraryscores arbitscores;
  unsigned int minlengthTSD,
               maxlengthTSD,
               allowedmismatches;
  unsigned int vicinity;
  bool bestoverlaps,
       nooverlaps,
       fastaoutput,
       fastaoutputinnerregion,
       gff3output,
       longoutput,
       scan,
       verbosemode;
  GtOption *optionmotif,
           *optionmotifmis,
           *optionoverlaps,
           *optionout,
           *optionoutinner,
           *optiongff3;
  GtLTRFourCharMotif *motif;
} LTRharvestArguments;

static void gt_ltrharvest_showopts(const LTRharvestArguments *lo)
{
  printf("# user defined options and values:\n");
  if (lo->verbosemode)
  {
    printf("#   verbosemode: On\n");
  }
  else
  {
    printf("#   verbosemode: Off\n");
  }
  printf("#   indexname: %s\n", gt_str_get(lo->str_indexname));
  if (lo->fastaoutput)
  {
    printf("#   outputfile: %s\n", gt_str_get(lo->str_fastaoutputfilename));
  }
  if (lo->fastaoutputinnerregion)
  {
    printf("#   outputfile inner region: %s\n",
        gt_str_get(lo->str_fastaoutputfilenameinnerregion));
  }
  if (lo->gff3output)
  {
    printf("#   outputfile gff3 format: %s\n",
           gt_str_get(lo->str_gff3filename));
  }
  printf("#   xdropbelowscore: %d\n", lo->xdropbelowscore);
  printf("#   similaritythreshold: %.2f\n", lo->similaritythreshold);
  printf("#   minseedlength: %lu\n", lo->minseedlength);
  printf("#   matchscore: %d\n", lo->arbitscores.mat);
  printf("#   mismatchscore: %d\n", lo->arbitscores.mis);
  printf("#   insertionscore: %d\n", lo->arbitscores.ins);
  printf("#   deletionscore: %d\n", lo->arbitscores.del);
  printf("#   minLTRlength: %lu\n",  lo->minltrlength);
  printf("#   maxLTRlength: %lu\n",  lo->maxltrlength);
  printf("#   minLTRdistance: %lu\n",  lo->mindistance);
  printf("#   maxLTRdistance: %lu\n",  lo->maxdistance);
  if (lo->nooverlaps)
  {
    printf("#   overlaps: no\n");
  }
  else
  {
    if (lo->bestoverlaps)
    {
      printf("#   overlaps: best\n");
    }
    else
    {
      printf("#   overlaps: all\n");
    }
  }
  printf("#   minTSDlength: %u\n",  lo->minlengthTSD);
  printf("#   maxTSDlength: %u\n",  lo->maxlengthTSD);
  printf("#   palindromic motif: %s\n", gt_str_get(lo->str_motif));
  printf("#   motifmismatchesallowed: %u\n", lo->allowedmismatches);
  printf("#   vicinity: %u nt\n", lo->vicinity);
  if (lo->searchrange.start != 0 ||
      lo->searchrange.end != 0)
  {
    printf("# ltrsearchseqrange=(%lu,%lu)\n",
          lo->searchrange.start,
          lo->searchrange.end);
  }
}

static void* gt_ltrharvest_arguments_new(void)
{
  LTRharvestArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->ofi = gt_output_file_info_new();
  arguments->str_motif = gt_str_new();
  arguments->str_overlaps = gt_str_new();
  arguments->str_indexname = gt_str_new();
  arguments->str_fastaoutputfilename = gt_str_new();
  arguments->str_fastaoutputfilenameinnerregion = gt_str_new();
  arguments->str_gff3filename = gt_str_new();
  arguments->motif = gt_calloc((size_t) 1, sizeof (GtLTRFourCharMotif));
  return arguments;
}

static void gt_ltrharvest_arguments_delete(void *tool_arguments)
{
  LTRharvestArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_output_file_info_delete(arguments->ofi);
  gt_str_delete(arguments->str_motif);
  gt_str_delete(arguments->str_overlaps);
  gt_str_delete(arguments->str_indexname);
  gt_str_delete(arguments->str_fastaoutputfilename);
  gt_str_delete(arguments->str_fastaoutputfilenameinnerregion);
  gt_str_delete(arguments->str_gff3filename);
  gt_option_delete(arguments->optiongff3);
  gt_option_delete(arguments->optionoutinner);
  gt_option_delete(arguments->optionout);
  gt_option_delete(arguments->optionmotif);
  gt_option_delete(arguments->optionmotifmis);
  gt_option_delete(arguments->optionoverlaps);
  gt_free(arguments->motif);
  gt_free(arguments);
}

static GtOptionParser* gt_ltrharvest_option_parser_new(void *tool_arguments)
{
  LTRharvestArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *optionindex,
           *optionltrsearchseqrange,
           *optionseed,
           *optionminlenltr,
           *optionmaxlenltr,
           *optionmindistltr,
           *optionmaxdistltr,
           *optionmintsd,
           *optionmaxtsd,
           *optionsimilar,
           *optionmotif,
           *optionmotifmis,
           *optionvic,
           *optionoverlaps,
           *optionxdrop,
           *optionmat,
           *optionmis,
           *optionins,
           *optiondel,
           *optionv,
           *optionoffset,
           *optionlongoutput,
           *optionout,
           *optionoutinner,
           *optiongff3,
           *optionscan;
  GtRange default_ltrsearchseqrange = {0,0};
  static const char *overlaps[] = {
    "best", /* the default */
    "no",
    "all",
    NULL
  };

  gt_assert(arguments);

  op = gt_option_parser_new("[option ...] -index <indexname>",
                            "Predict LTR retrotransposons.");

  /* -index */
  optionindex = gt_option_new_string("index",
                             "specify the name of the enhanced suffix "
                             "array index (mandatory)",
                             arguments->str_indexname, NULL);
  gt_option_is_mandatory(optionindex);
  gt_option_parser_add_option(op, optionindex);

  /* -range */
  optionltrsearchseqrange
    = gt_option_new_range("range",
                          "specify range in the input sequence(s) in which LTR "
                          "pairs are searched",
                          &arguments->searchrange,
                          &default_ltrsearchseqrange);
  gt_option_parser_add_option(op, optionltrsearchseqrange);

  /* -seed */
  optionseed = gt_option_new_ulong_min("seed",
                               "specify minimum seed length for"
                               " exact repeats",
                               &arguments->minseedlength,
                               30UL,
                               1UL);
  gt_option_parser_add_option(op, optionseed);

  /* -minlenltr */
  optionminlenltr = gt_option_new_ulong_min_max("minlenltr",
                               "specify minimum length for each LTR",
                               &arguments->minltrlength,
                               100UL,
                               1UL,
                               GT_UNDEF_ULONG);
  gt_option_parser_add_option(op, optionminlenltr);

  /* -maxlenltr */
  optionmaxlenltr = gt_option_new_ulong_min_max("maxlenltr",
                               "specify maximum length for each LTR",
                               &arguments->maxltrlength,
                               1000UL,
                               1UL,
                               GT_UNDEF_ULONG);
  gt_option_parser_add_option(op, optionmaxlenltr);

  /* -mindistltr */
  optionmindistltr = gt_option_new_ulong_min_max("mindistltr",
                               "specify minimum distance of "
                               "LTR startpositions",
                               &arguments->mindistance,
                               1000UL,
                               1UL,
                               GT_UNDEF_ULONG);
  gt_option_parser_add_option(op, optionmindistltr);

  /* -maxdistltr */
  optionmaxdistltr = gt_option_new_ulong_min_max("maxdistltr",
                               "specify maximum distance of "
                               "LTR startpositions",
                               &arguments->maxdistance,
                               15000UL,
                               1UL,
                               GT_UNDEF_ULONG);
  gt_option_parser_add_option(op, optionmaxdistltr);

  /* -similar */
  optionsimilar = gt_option_new_double_min_max("similar",
                               "specify similaritythreshold in "
                               "range [1..100%]",
                               &arguments->similaritythreshold,
                               (double) 85.0,
                               (double) 0.0,
                               100.0);
  gt_option_parser_add_option(op, optionsimilar);

  /* -mintsd */
  optionmintsd = gt_option_new_uint_min_max("mintsd",
                              "specify minimum length for each TSD",
                               &arguments->minlengthTSD,
                               4U,
                               0,
                               GT_UNDEF_UINT);
  gt_option_parser_add_option(op, optionmintsd);

  /* -maxtsd */
  optionmaxtsd = gt_option_new_uint_min_max("maxtsd",
                              "specify maximum length for each TSD",
                               &arguments->maxlengthTSD,
                               20U,
                               0,
                               GT_UNDEF_UINT);
  gt_option_parser_add_option(op, optionmaxtsd);

  optionmotif = gt_option_new_string("motif",
                             "specify 2 nucleotides startmotif + "
                             "2 nucleotides endmotif: ****",
                             arguments->str_motif, NULL);
  gt_option_parser_add_option(op, optionmotif);
  arguments->optionmotif = gt_option_ref(optionmotif);

  /* -motifmis */
  optionmotifmis = gt_option_new_uint_min_max("motifmis",
                             "specify maximum number of "
                             "mismatches in motif [0,3]",
                             &arguments->allowedmismatches,
                             4U,
                             0,
                             3U);
  gt_option_parser_add_option(op, optionmotifmis);
  arguments->optionmotifmis = gt_option_ref(optionmotifmis);

  /* -vic */
  optionvic = gt_option_new_uint_min_max("vic",
                        "specify the number of nucleotides (to the left and "
                        "to the right) that will be searched "
                        "for TSDs and/or motifs around 5' and 3' boundary "
                        "of predicted LTR retrotransposons",
                        &arguments->vicinity,
                        60U,
                        1U,
                        500U);
  gt_option_parser_add_option(op, optionvic);

  /* -overlaps */
  optionoverlaps = gt_option_new_choice("overlaps",
               "specify no|best|all",
               arguments->str_overlaps,
               overlaps[0],
               overlaps);
  gt_option_parser_add_option(op, optionoverlaps);
  arguments->optionoverlaps = gt_option_ref(optionoverlaps);

  /* -xdrop */
  optionxdrop = gt_option_new_int_min("xdrop",
                        "specify xdropbelowscore for extension-alignment",
                        &arguments->xdropbelowscore,
                        5,
                        0);
  gt_option_parser_add_option(op, optionxdrop);

  /* -mat */
  optionmat = gt_option_new_int_min("mat",
                        "specify matchscore for extension-alignment",
                        &arguments->arbitscores.mat,
                        2,
                        1);
  gt_option_parser_add_option(op, optionmat);

  /* -mis */
  optionmis = gt_option_new_int_max("mis",
                        "specify mismatchscore for extension-alignment",
                        &arguments->arbitscores.mis,
                        -2,
                        -1);
  gt_option_parser_add_option(op, optionmis);

  /* -ins */
  optionins = gt_option_new_int_max("ins",
                        "specify insertionscore for extension-alignment",
                        &arguments->arbitscores.ins,
                        -3,
                        -1);
  gt_option_parser_add_option(op, optionins);

  /* -del */
  optiondel = gt_option_new_int_max("del",
                        "specify deletionscore for extension-alignment",
                        &arguments->arbitscores.del,
                        -3,
                        -1);
  gt_option_parser_add_option(op, optiondel);

  /* -v */
  optionv = gt_option_new_bool("v",
                           "verbose mode",
                           &arguments->verbosemode,
                           false);
  gt_option_parser_add_option(op, optionv);

  /* -longoutput */
  optionlongoutput = gt_option_new_bool("longoutput",
                           "additional motif/TSD output",
                           &arguments->longoutput,
                           false);
  gt_option_parser_add_option(op, optionlongoutput);

  /* -out */
  arguments->fastaoutput = false;
  optionout = gt_option_new_string("out",
                             "specify FASTA outputfilename",
                             arguments->str_fastaoutputfilename, NULL);
  gt_option_parser_add_option(op, optionout);
  arguments->optionout = gt_option_ref(optionout);

  /* -outinner */
  arguments->fastaoutputinnerregion = false;
  optionoutinner = gt_option_new_string("outinner",
                             "specify FASTA outputfilename for inner regions",
                             arguments->str_fastaoutputfilenameinnerregion,
                              NULL);
  gt_option_parser_add_option(op, optionoutinner);
  arguments->optionoutinner = gt_option_ref(optionoutinner);

  /* -gff3 */
  arguments->gff3output = false;
  optiongff3 = gt_option_new_string("gff3",
                                    "specify GFF3 outputfilename",
                                    arguments->str_gff3filename, NULL);
  gt_option_parser_add_option(op, optiongff3);
  arguments->optiongff3 = gt_option_ref(optiongff3);

  /* -offset */
  optionoffset = gt_option_new_ulong("offset",
                                     "offset added to GFF3 coordinates",
                                     &arguments->offset,
                                     0UL);
  gt_option_parser_add_option(op, optionoffset);
  gt_option_is_extended_option(optionoffset);

  /* -scan */
  optionscan = gt_option_new_bool("scan",
                                  "scan the index sequentially instead of "
                                  "mapping it into memory entirely",
                                  &arguments->scan,
                                  true);
  gt_option_parser_add_option(op, optionscan);
  gt_option_is_extended_option(optionscan);

  /* implications */
  gt_option_imply(optionmaxtsd, optionmintsd);
  gt_option_imply(optionmotifmis, optionmotif);

  gt_option_imply_either_2(optionlongoutput, optionmintsd, optionmotif);

  gt_option_parser_refer_to_manual(op);

  /* gt_output_file_register_options(op, &arguments->outfp, arguments->ofi); */

  return op;
}

static int gt_ltrharvest_arguments_check(GT_UNUSED int rest_argc,
                                           void *tool_arguments,
                                           GtError *err)
{
  LTRharvestArguments *arguments = tool_arguments;
  int had_err = 0;

  /* init */
  arguments->motif->firstleft   = (GtUchar) 't';
  arguments->motif->secondleft  = (GtUchar) 'g';
  arguments->motif->firstright  = (GtUchar) 'c';
  arguments->motif->secondright = (GtUchar) 'a';
  arguments->motif->allowedmismatches = arguments->allowedmismatches;

  if (gt_option_is_set(arguments->optionmotif))
  {
    if (gt_str_length(arguments->str_motif) != 4UL)
    {
      gt_error_set(err,
          "argument of -motif has not exactly 4 characters");
      had_err = -1;
    }
    if (!had_err) {
      arguments->motif->firstleft =
        (GtUchar) gt_str_get(arguments->str_motif)[0];
      arguments->motif->secondleft =
        (GtUchar) gt_str_get(arguments->str_motif)[1];
      arguments->motif->firstright =
        (GtUchar)  gt_str_get(arguments->str_motif)[2];
      arguments->motif->secondright =
        (GtUchar)  gt_str_get(arguments->str_motif)[3];
      /* default if motif specified */
      if (!gt_option_is_set(arguments->optionmotifmis))
      {
        arguments->motif->allowedmismatches = 0;
        arguments->allowedmismatches = 0;
      }
    }
  }
  if (!had_err) {
    if (gt_option_is_set(arguments->optionoverlaps))
    {
      if (strcmp(gt_str_get(arguments->str_overlaps), "no") == 0)
      {
        arguments->bestoverlaps = false;
        arguments->nooverlaps = true;
      }
      else if (strcmp(gt_str_get(arguments->str_overlaps), "best") == 0 )
      {
        arguments->bestoverlaps = true;
        arguments->nooverlaps = false;
      }
      else if (strcmp(gt_str_get(arguments->str_overlaps), "all") == 0 )
      {
        arguments->bestoverlaps = false;
        arguments->nooverlaps = false;
      }
      else
      {
        gt_assert(0); /* cannot happen */
      }
    }
    else
    {
      /* default is "best" */
      arguments->bestoverlaps = true;  /* take best prediction
                                          if overlap occurs, default */
      arguments->nooverlaps = false; /* overlapping predictions (not) allowed*/
    }
    /* if FASTA output is set */
    if (gt_option_is_set(arguments->optionout))
    {
      arguments->fastaoutput = true;
    }
    /* if FASTA output inner region is set */
    if (gt_option_is_set(arguments->optionoutinner))
    {
      arguments->fastaoutputinnerregion = true;
    }
    /* if GFF3 output is set */
    if (gt_option_is_set(arguments->optiongff3))
    {
      arguments->gff3output = true;
    }
  }
  return had_err;
}

static void gt_ltrharvest_showargsline(int argc, const char **argv)
{
  int i;
  gt_assert(argv && argc >= 1);
  printf("# args=");
  for (i=1; i<argc; i++) {
    printf("%s", argv[i]);
    if (i != argc-1) printf(" ");
  }
  printf("\n");
}

static int gt_ltrharvest_runner(GT_UNUSED int argc,
                                  GT_UNUSED const char **argv,
                                  GT_UNUSED int parsed_args,
                                  void *tool_arguments, GtError *err)
{
  GtNodeStream *ltrh_stream,
               *gff3_out_stream = NULL,
               *tabout_stream = NULL,
               *fasta_out_stream = NULL,
               *fasta_inner_out_stream = NULL,
               *last_stream;
  GtNodeVisitor *tabout_visitor = NULL;
  GtFile *gff3file = NULL,
         *fastaoutfile = NULL,
         *fastainneroutfile = NULL;
  LTRharvestArguments *arguments = tool_arguments;
  const GtEncseq *encseq = NULL;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  /* create a LTRharvest stream */
  ltrh_stream = gt_ltrharvest_stream_new(arguments->str_indexname,
                                         arguments->searchrange,
                                         arguments->minseedlength,
                                         arguments->minltrlength,
                                         arguments->maxltrlength,
                                         arguments->mindistance,
                                         arguments->maxdistance,
                                         arguments->similaritythreshold,
                                         arguments->xdropbelowscore,
                                         arguments->arbitscores,
                                         arguments->motif,
                                         arguments->verbosemode,
                                         arguments->nooverlaps,
                                         arguments->bestoverlaps,
                                         arguments->scan,
                                         arguments->offset,
                                         arguments->minlengthTSD,
                                         arguments->maxlengthTSD,
                                         (unsigned long) arguments->vicinity,
                                         err);
  if (ltrh_stream == NULL)
    return -1;
  last_stream = ltrh_stream;

  encseq = gt_ltrharvest_stream_get_encseq(ltrh_stream);

  /* set visitors according to output parameter */
  if (!had_err && arguments->longoutput) {
    tabout_visitor = gt_ltrharvest_tabout_visitor_new_longoutput(encseq);
  } else {
    tabout_visitor = gt_ltrharvest_tabout_visitor_new();
  }

  /* create tabular output stream (traditional LTRharvest format) */
  tabout_stream = gt_ltrharvest_tabout_stream_new(last_stream, tabout_visitor);
  last_stream = tabout_stream;

  /* attach GFF3 output stream if requested */
  if (!had_err && arguments->gff3output) {
    gff3file = gt_file_open(GT_FILE_MODE_UNCOMPRESSED,
                            gt_str_get(arguments->str_gff3filename),
                            "w+",
                            err);
    if (gff3file == NULL) {
      had_err = -1;
    } else {
      gff3_out_stream = gt_gff3_out_stream_new(last_stream, gff3file);
      last_stream = gff3_out_stream;
    }
  }

  /* attach FASTA output stream if requested */
  if (!had_err && arguments->fastaoutput) {
    fastaoutfile = gt_file_open(GT_FILE_MODE_UNCOMPRESSED,
                                gt_str_get(arguments->str_fastaoutputfilename),
                                "w+",
                                err);
    if (fastaoutfile == NULL) {
      had_err = -1;
    } else {
      fasta_out_stream = gt_ltrharvest_fasta_out_stream_new(last_stream,
                                                            false,
                                                            encseq,
                                                            60UL,
                                                            fastaoutfile);
      last_stream = fasta_out_stream;
    }
  }

  /* attach FASTA inner region output stream if requested */
  if (!had_err && arguments->fastaoutputinnerregion) {
    fastainneroutfile = gt_file_open(GT_FILE_MODE_UNCOMPRESSED,
                      gt_str_get(arguments->str_fastaoutputfilenameinnerregion),
                      "w+", err);
    if (fastainneroutfile == NULL) {
      had_err = -1;
    } else {
      fasta_inner_out_stream = gt_ltrharvest_fasta_out_stream_new(last_stream,
                                                             true,
                                                             encseq,
                                                             60UL,
                                                             fastainneroutfile);
      last_stream = fasta_inner_out_stream;
    }
  }

  /* output arguments line */
  gt_ltrharvest_showargsline(argc, argv);

  /* show long parameters */
  if (arguments->verbosemode)
    gt_ltrharvest_showopts(arguments);

  /* print tabular output header */
  if (!had_err) {
    if (arguments->longoutput) {
      gt_ltrharvest_tabout_stream_printlongheader(arguments->minlengthTSD > 1U,
                                             arguments->allowedmismatches < 4U);
    } else {
      gt_ltrharvest_tabout_stream_printshortheader();
    }
  }

  /* pull the features through the stream and free them afterwards */
  if (!had_err)
    had_err = gt_node_stream_pull(last_stream, err);

  /* free */
  gt_node_stream_delete(gff3_out_stream);
  gt_node_stream_delete(tabout_stream);
  gt_node_stream_delete(fasta_out_stream);
  gt_node_stream_delete(fasta_inner_out_stream);
  gt_node_stream_delete(ltrh_stream);

  if (gff3file != NULL) gt_file_delete(gff3file);
  if (fastaoutfile != NULL) gt_file_delete(fastaoutfile);
  if (fastainneroutfile != NULL) gt_file_delete(fastainneroutfile);

  return had_err;
}

GtTool* gt_ltrharvest(void)
{
  return gt_tool_new(gt_ltrharvest_arguments_new,
                     gt_ltrharvest_arguments_delete,
                     gt_ltrharvest_option_parser_new,
                     gt_ltrharvest_arguments_check,
                     gt_ltrharvest_runner);
}
