/*
  Copyright (c) 2007 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Malte Mader <mmader@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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
#include <cairo.h>
#include "core/cstr.h"
#include "core/fileutils.h"
#include "core/gtdatapath.h"
#include "core/option.h"
#include "core/splitter.h"
#include "core/undef.h"
#include "core/versionfunc.h"
#include "extended/add_introns_stream.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream.h"
#include "annotationsketch/canvas.h"
#include "annotationsketch/diagram.h"
#include "annotationsketch/feature_index.h"
#include "annotationsketch/feature_stream.h"
#include "annotationsketch/gt_sketch.h"
#include "annotationsketch/image_info.h"
#include "annotationsketch/style.h"

typedef struct {
  bool pipe,
       verbose,
       addintrons,
       showrecmaps;
  Str *seqid, *format;
  unsigned long start,
                end;
  unsigned int width;
} AnnotationSketchArguments;

static OPrval parse_options(int *parsed_args,
                            AnnotationSketchArguments *arguments,
                            int argc, const char **argv, GT_Error *err)
{
  OptionParser *op;
  Option  *option, *option2;
  OPrval oprval;
  bool force;
  static const char *formats[] = { "png",
#ifdef CAIRO_HAS_PDF_SURFACE
  "pdf",
#endif
#ifdef CAIRO_HAS_SVG_SURFACE
  "svg",
#endif
#ifdef CAIRO_HAS_PS_SURFACE
  "ps",
#endif
  NULL };
  gt_error_check(err);

  /* init */
  op = option_parser_new("[option ...] image_file [GFF3_file ...]",
                         "Create graphical representations of "
                         "GFF3 annotation files.");

  /* -v */
  option = option_new_verbose(&arguments->verbose);
  option_parser_add_option(op, option);

  /* -pipe */
  option = option_new_bool("pipe", "use pipe mode (i.e., show all gff3 "
                           "features on stdout)", &arguments->pipe, false);
  option_parser_add_option(op, option);

  /* -force */
  option = option_new_bool("force", "force writing to output file", &force,
                           false);
  option_parser_add_option(op, option);

  /* -seqid */
  option = option_new_string("seqid", "sequence region identifier\n"
                                      "default: first one in file",
                            arguments->seqid, NULL);
  option_parser_add_option(op, option);
  option_hide_default(option);

  /* -start */
  option = option_new_ulong_min("start", "start position\n"
                                         "default: first region start",
                            &arguments->start, UNDEF_ULONG, 1);
  option_parser_add_option(op, option);
  option_hide_default(option);

  /* -end */
  option2 = option_new_ulong("end", "end position\ndefault: last region end",
                            &arguments->end, UNDEF_ULONG);
  option_parser_add_option(op, option2);
  /* -start and -end must be given together */
  option_imply(option, option2);
  option_imply(option2, option);
  option_hide_default(option2);

  /* -width */
  option = option_new_uint_min("width", "target image width", &arguments->width,
                               800, 1);
  option_parser_add_option(op, option);

  /* -format */
  option = option_new_choice("format", "output graphics format\n"
                                       "choose from png"
#ifdef CAIRO_HAS_PDF_SURFACE
                                       "|pdf"
#endif
#ifdef CAIRO_HAS_SVG_SURFACE
                                       "|svg"
#endif
#ifdef CAIRO_HAS_PS_SURFACE
                                       "|ps"
#endif
                                       "",
                             arguments->format, formats[0], formats);
  option_parser_add_option(op, option);

  /* -addintrons */
  option = option_new_bool("addintrons", "add intron features between "
                           "existing exon features (before drawing)",
                           &arguments->addintrons, false);
  option_parser_add_option(op, option);

  /* -showrecmaps */
  option = option_new_bool("showrecmaps", "show RecMaps after image creation",
                           &arguments->showrecmaps, false);
  option_is_development_option(option);
  option_parser_add_option(op, option);

  /* set contact mailaddress */
  option_parser_set_mailaddress(op, "<ssteinbiss@stud.zbh.uni-hamburg.de>");

  /* parse options */
  option_parser_set_min_args(op, 1);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);

  if (oprval == OPTIONPARSER_OK && !force && file_exists(argv[*parsed_args])) {
    gt_error_set(err, "file \"%s\" exists already. use option -force to "
                   "overwrite", argv[*parsed_args]);
    oprval = OPTIONPARSER_ERROR;
  }

  /* free */
  option_parser_delete(op);

  return oprval;
}

int gt_sketch(int argc, const char **argv, GT_Error *err)
{
  GenomeStream *gff3_in_stream = NULL,
               *add_introns_stream = NULL,
               *gff3_out_stream = NULL,
               *feature_stream = NULL,
               *last_stream;
  AnnotationSketchArguments arguments;
  GT_GenomeNode *gn = NULL;
  GT_FeatureIndex *features = NULL;
  int parsed_args, had_err=0;
  const char *file, *seqid = NULL;
  GT_Range qry_range, sequence_region_range;
  GT_Array *results = NULL;
  GT_Style *sty = NULL;
  Str *gt_style_file = NULL;
  Str *prog;
  GT_Diagram *d = NULL;
  GT_ImageInfo* ii = NULL;
  GT_Canvas *canvas = NULL;

  gt_error_check(err);

  /* option parsing */
  arguments.seqid = str_new();
  arguments.format = str_new();
  switch (parse_options(&parsed_args, &arguments, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
      str_delete(arguments.seqid);
      str_delete(arguments.format);
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      str_delete(arguments.seqid);
      str_delete(arguments.format);
      return 0;
  }

  /* save name of output file */
  file = argv[parsed_args];

  /* check for correct order: range end < range start */
  if (!had_err &&
      arguments.start != UNDEF_ULONG &&
      arguments.end != UNDEF_ULONG &&
      !(arguments.start < arguments.end)) {
    gt_error_set(err, "start of query range (%lu) must be before "
                   "end of query range (%lu)",
              arguments.start, arguments.end);
    had_err = -1;
  }

  if (!had_err) {
    /* create feature index */
    features = gt_feature_index_new();
    parsed_args++;

    /* create a gff3 input stream */
    gff3_in_stream = gff3_in_stream_new_unsorted(argc - parsed_args,
                                                 argv + parsed_args,
                                                 arguments.verbose, false);
    last_stream = gff3_in_stream;

    /* create add introns stream if -addintrons was used */
    if (arguments.addintrons) {
      add_introns_stream = add_introns_stream_new(last_stream);
      last_stream = add_introns_stream;
    }

    /* create gff3 output stream if -pipe was used */
    if (arguments.pipe) {
      gff3_out_stream = gff3_out_stream_new(last_stream, NULL);
      last_stream = gff3_out_stream;
    }

    /* create feature stream */
    feature_stream = feature_stream_new(last_stream, features);

    /* pull the features through the stream and free them afterwards */
    while (!(had_err = genome_stream_next_tree(feature_stream, &gn, err)) &&
           gn) {
      genome_node_rec_delete(gn);
    }

    genome_stream_delete(feature_stream);
    genome_stream_delete(gff3_out_stream);
    genome_stream_delete(add_introns_stream);
    genome_stream_delete(gff3_in_stream);
  }

  /* if seqid is empty, take first one added to index */
  if (!had_err && strcmp(str_get(arguments.seqid),"") == 0) {
    seqid = gt_feature_index_get_first_seqid(features);
    if (seqid == NULL) {
      gt_error_set(err, "GFF input file must contain a sequence region!");
      had_err = -1;
    }
  }
  else if (!had_err && !gt_feature_index_has_seqid(features,
                                                str_get(arguments.seqid))) {
    gt_error_set(err, "sequence region '%s' does not exist in GFF input file",
              str_get(arguments.seqid));
    had_err = -1;
  }
  else if (!had_err)
    seqid = str_get(arguments.seqid);

  results = gt_array_new(sizeof (GT_GenomeNode*));
  if (!had_err) {
    gt_feature_index_get_range_for_seqid(features, &sequence_region_range,
                                         seqid);
    qry_range.start = (arguments.start == UNDEF_ULONG ?
                         sequence_region_range.start :
                         arguments.start);
    qry_range.end   = (arguments.end == UNDEF_ULONG ?
                         sequence_region_range.end :
                         arguments.end);
  }

  if (!had_err) {
    if (arguments.verbose)
      fprintf(stderr, "# of results: %lu\n", gt_array_size(results));

    /* find and load style file */
    prog = str_new();
    str_append_cstr_nt(prog, argv[0], cstr_length_up_to_char(argv[0], ' '));
    gt_style_file = gtdata_get_path(str_get(prog), err);
    str_delete(prog);
    str_append_cstr(gt_style_file, "/sketch/default.style");
    if (!(sty = gt_style_new(arguments.verbose, err)))
      had_err = -1;
    if (!had_err && file_exists(str_get(gt_style_file)))
      had_err = gt_style_load_file(sty, str_get(gt_style_file), err);
  }

  if (!had_err) {
    /* create and write image file */
    d = gt_diagram_new(features, seqid, &qry_range, sty);
    ii = gt_image_info_new();
    if (strcmp(str_get(arguments.format),"pdf")==0)
      canvas = gt_canvas_new(sty, GRAPHICS_PDF, arguments.width, ii);
    else if (strcmp(str_get(arguments.format),"ps")==0)
      canvas = gt_canvas_new(sty, GRAPHICS_PS, arguments.width, ii);
    else if (strcmp(str_get(arguments.format),"svg")==0)
      canvas = gt_canvas_new(sty, GRAPHICS_SVG, arguments.width, ii);
    else
      canvas = gt_canvas_new(sty, GRAPHICS_PNG, arguments.width, ii);
    gt_diagram_sketch(d, canvas);
    if (arguments.showrecmaps) {
      unsigned long i;
      const GT_RecMap *rm;
      for (i = 0; i < gt_image_info_num_of_recmaps(ii) ;i++) {
        GenomeFeatureType *type;
        char buf[BUFSIZ];
        rm = gt_image_info_get_recmap(ii, i);
        gt_recmap_format_html_imagemap_coords(rm, buf, BUFSIZ);
        type = genome_feature_get_type(rm->gf);
        printf("%s, %s\n", buf, genome_feature_type_get_cstr(type));
      }
    }
    had_err = gt_canvas_to_file(canvas, file, err);
  }

  /* free */
  gt_canvas_delete(canvas);
  gt_image_info_delete(ii);
  gt_style_delete(sty);
  str_delete(gt_style_file);
  gt_diagram_delete(d);
  str_delete(arguments.seqid);
  str_delete(arguments.format);
  gt_array_delete(results);
  gt_feature_index_delete(features);

  return had_err;
}
