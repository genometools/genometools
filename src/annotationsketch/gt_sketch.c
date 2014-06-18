/*
  Copyright (c) 2007-2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007      Malte Mader <mader@zbh.uni-hamburg.de>
  Copyright (c) 2007      Christin Schaerfer <schaerfer@zbh.uni-hamburg.de>
  Copyright (c) 2007-2013 Center for Bioinformatics, University of Hamburg

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

#include <cairo.h>
#include <string.h>
#include "core/cstr_api.h"
#include "core/fileutils_api.h"
#include "core/gtdatapath.h"
#include "core/option_api.h"
#include "core/output_file_api.h"
#include "core/ma.h"
#include "core/splitter.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/versionfunc.h"
#include "core/warning_api.h"
#include "extended/add_introns_stream_api.h"
#include "extended/bed_in_stream.h"
#include "extended/feature_index_memory_api.h"
#include "extended/feature_stream_api.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream_api.h"
#include "extended/gtf_in_stream.h"
#include "extended/sort_stream_api.h"
#include "annotationsketch/block.h"
#include "annotationsketch/canvas_api.h"
#include "annotationsketch/canvas_cairo_file.h"
#include "annotationsketch/diagram.h"
#include "annotationsketch/gt_sketch.h"
#include "annotationsketch/image_info.h"
#include "annotationsketch/layout.h"
#include "annotationsketch/style.h"

typedef struct {
  bool pipe,
       verbose,
       addintrons,
       showrecmaps,
       flattenfiles,
       unsafe,
       force,
       use_streams;
  GtStr *seqid, *format, *stylefile, *input;
  GtUword start,
                end;
  unsigned int width;
} GtSketchArguments;

static void* gt_sketch_arguments_new(void)
{
  GtSketchArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->seqid = gt_str_new();
  arguments->format = gt_str_new();
  arguments->input = gt_str_new();
  arguments->stylefile = gt_str_new();
  return arguments;
}

static void gt_sketch_arguments_delete(void *tool_arguments)
{
  GtSketchArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->seqid);
  gt_str_delete(arguments->format);
  gt_str_delete(arguments->input);
  gt_str_delete(arguments->stylefile);
  gt_free(arguments);
}

static GtOptionParser* gt_sketch_option_parser_new(void *tool_arguments)
{
  GtSketchArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *option2;
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
    NULL
  };
  static const char *inputs[] = {
    "gff",
    "bed",
    "gtf",
    NULL
  };
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] image_file [GFF3_file ...]",
                            "Create graphical representation of GFF3 "
                            "annotation files.");

  /* -pipe */
  option = gt_option_new_bool("pipe", "use pipe mode (i.e., show all gff3 "
                              "features on stdout)", &arguments->pipe, false);
  gt_option_parser_add_option(op, option);

  /* -flattenfiles */
  option = gt_option_new_bool("flattenfiles", "do not group tracks by source "
                              "file name and remove file names from track "
                              "description", &arguments->flattenfiles, false);
  gt_option_parser_add_option(op, option);

  /* -seqid */
  option = gt_option_new_string("seqid", "sequence region identifier\n"
                                      "default: first one in file",
                            arguments->seqid, NULL);
  gt_option_parser_add_option(op, option);
  gt_option_hide_default(option);

  /* -start */
  option = gt_option_new_uword_min("start", "start position\n"
                                         "default: first region start",
                            &arguments->start, GT_UNDEF_UWORD, 1);
  gt_option_parser_add_option(op, option);
  gt_option_hide_default(option);

  /* -end */
  option2 = gt_option_new_uword("end", "end position\ndefault: last region end",
                            &arguments->end, GT_UNDEF_UWORD);
  gt_option_parser_add_option(op, option2);
  /* -start and -end must be given together */
  gt_option_imply(option, option2);
  gt_option_imply(option2, option);
  gt_option_hide_default(option2);

  /* -width */
  option = gt_option_new_uint_min("width", "target image width (in pixel)",
                                  &arguments->width,
                                  800, 1);
  gt_option_parser_add_option(op, option);

  /* -style */
  option = gt_option_new_string("style", "style file to use",
                                arguments->stylefile,
                                gt_str_get(arguments->stylefile));
  gt_option_parser_add_option(op, option);

  /* -format */
  option = gt_option_new_choice("format", "output graphics format\n"
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
  gt_option_parser_add_option(op, option);

  /* -input */
  option = gt_option_new_choice("input", "input data format\n"
                                       "choose from gff|bed|gtf",
                             arguments->input, inputs[0], inputs);
  gt_option_parser_add_option(op, option);

  /* -addintrons */
  option = gt_option_new_bool("addintrons", "add intron features between "
                              "existing exon features (before drawing)",
                              &arguments->addintrons, false);
  gt_option_parser_add_option(op, option);

    /* -unsafe */
  option = gt_option_new_bool("unsafe", "enable unsafe mode for style file",
                              &arguments->unsafe, false);
  gt_option_parser_add_option(op, option);

  /* -showrecmaps */
  option = gt_option_new_bool("showrecmaps",
                              "show RecMaps after image creation",
                              &arguments->showrecmaps, false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -streams */
  option = gt_option_new_bool("streams", "use streams to write data to file",
                              &arguments->use_streams, false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  /* -force */
  option = gt_option_new_bool(GT_FORCE_OPT_CSTR, "force writing to output file",
                              &arguments->force, false);
  gt_option_parser_add_option(op, option);

  gt_option_parser_set_min_args(op, 1);

  return op;
}

static int gt_sketch_arguments_check(GT_UNUSED int rest_argc,
                                     void *tool_arguments,
                                     GT_UNUSED GtError *err)
{
  GtSketchArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  if (arguments->start != GT_UNDEF_UWORD &&
      arguments->end != GT_UNDEF_UWORD &&
      !(arguments->start < arguments->end)) {
    gt_error_set(err, "start of query range ("GT_WU") must be before "
                      "end of query range ("GT_WU")",
                      arguments->start, arguments->end);
    had_err = -1;
  }

  return had_err;
}

/* this track selector function is used to disregard file names in track
   identifiers */
static void flattened_file_track_selector(GtBlock *block, GtStr *result,
                                   GT_UNUSED void *data)
{
  gt_assert(block && result);
  gt_str_reset(result);
  gt_str_append_cstr(result, gt_block_get_type(block));
}

static int gt_sketch_runner(int argc, const char **argv, int parsed_args,
                              void *tool_arguments, GT_UNUSED GtError *err)
{
  GtSketchArguments *arguments = tool_arguments;
  GtNodeStream *in_stream = NULL,
               *add_introns_stream = NULL,
               *gff3_out_stream = NULL,
               *feature_stream = NULL,
               *sort_stream = NULL,
               *last_stream;
  GtFeatureIndex *features = NULL;
  const char *file;
  char *seqid = NULL;
  GtRange qry_range, sequence_region_range;
  GtArray *results = NULL;
  GtStyle *sty = NULL;
  GtStr *prog, *defaultstylefile = NULL;
  GtDiagram *d = NULL;
  GtLayout *l = NULL;
  GtImageInfo* ii = NULL;
  GtCanvas *canvas = NULL;
  GtUword height;
  bool has_seqid;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  prog = gt_str_new();
  gt_str_append_cstr_nt(prog, argv[0],
                        gt_cstr_length_up_to_char(argv[0], ' '));
  defaultstylefile = gt_get_gtdata_path(gt_str_get(prog), err);
  gt_str_delete(prog);
  if (!defaultstylefile)
    had_err = -1;
  if (!had_err) {
    gt_str_append_cstr(defaultstylefile, "/sketch/default.style");
  }

  file = argv[parsed_args];
  if (!had_err) {
    /* create feature index */
    features = gt_feature_index_memory_new();
    parsed_args++;

    /* create an input stream */
    if (strcmp(gt_str_get(arguments->input), "gff") == 0)
    {
      in_stream = gt_gff3_in_stream_new_unsorted(argc - parsed_args,
                                                 argv + parsed_args);
      if (arguments->verbose)
        gt_gff3_in_stream_show_progress_bar((GtGFF3InStream*) in_stream);
    } else if (strcmp(gt_str_get(arguments->input), "bed") == 0)
    {
      if (argc - parsed_args == 0)
        in_stream = gt_bed_in_stream_new(NULL);
      else
        in_stream = gt_bed_in_stream_new(argv[parsed_args]);
    } else if (strcmp(gt_str_get(arguments->input), "gtf") == 0)
    {
      if (argc - parsed_args == 0)
        in_stream = gt_gtf_in_stream_new(NULL);
      else
        in_stream = gt_gtf_in_stream_new(argv[parsed_args]);
    }
    last_stream = in_stream;

    /* create add introns stream if -addintrons was used */
    if (arguments->addintrons) {
      sort_stream = gt_sort_stream_new(last_stream);
      add_introns_stream = gt_add_introns_stream_new(sort_stream);
      last_stream = add_introns_stream;
    }

    /* create gff3 output stream if -pipe was used */
    if (arguments->pipe) {
      gff3_out_stream = gt_gff3_out_stream_new(last_stream, NULL);
      last_stream = gff3_out_stream;
    }

    /* create feature stream */
    feature_stream = gt_feature_stream_new(last_stream, features);

    /* pull the features through the stream and free them afterwards */
    had_err = gt_node_stream_pull(feature_stream, err);

    gt_node_stream_delete(feature_stream);
    gt_node_stream_delete(gff3_out_stream);
    gt_node_stream_delete(sort_stream);
    gt_node_stream_delete(add_introns_stream);
    gt_node_stream_delete(in_stream);
  }

  if (!had_err) {
    had_err = gt_feature_index_has_seqid(features,
                                         &has_seqid,
                                         gt_str_get(arguments->seqid),
                                         err);
  }

  /* if seqid is empty, take first one added to index */
  if (!had_err && strcmp(gt_str_get(arguments->seqid),"") == 0) {
    seqid = gt_feature_index_get_first_seqid(features, err);
    if (seqid == NULL) {
      gt_error_set(err, "GFF input file must contain a sequence region!");
      had_err = -1;
    }
  }
  else if (!had_err && !has_seqid) {
    gt_error_set(err, "sequence region '%s' does not exist in GFF input file",
                 gt_str_get(arguments->seqid));
    had_err = -1;
  }
  else if (!had_err)
    seqid = gt_cstr_dup(gt_str_get(arguments->seqid));

  results = gt_array_new(sizeof (GtGenomeNode*));
  if (!had_err) {
    had_err = gt_feature_index_get_range_for_seqid(features,
                                                   &sequence_region_range,
                                                   seqid,
                                                   err);
  }
  if (!had_err) {
    qry_range.start = (arguments->start == GT_UNDEF_UWORD ?
                         sequence_region_range.start :
                         arguments->start);
    qry_range.end   = (arguments->end == GT_UNDEF_UWORD ?
                         sequence_region_range.end :
                         arguments->end);
  }

  if (!had_err) {
    if (arguments->verbose)
      fprintf(stderr, "# of results: "GT_WU"\n", gt_array_size(results));

    /* find and load style file */
    if (!(sty = gt_style_new(err)))
      had_err = -1;
    if (gt_str_length(arguments->stylefile) == 0) {
      gt_str_append_str(arguments->stylefile, defaultstylefile);
    } else {
      if (!had_err && gt_file_exists(gt_str_get(arguments->stylefile))) {
        if (arguments->unsafe)
          gt_style_unsafe_mode(sty);
      }
      else
      {
        had_err = -1;
        gt_error_set(err, "style file '%s' does not exist!",
                          gt_str_get(arguments->stylefile));
      }
    }
    if (!had_err)
      had_err = gt_style_load_file(sty, gt_str_get(arguments->stylefile), err);
  }

  if (!had_err) {
    /* create and write image file */
    if (!(d = gt_diagram_new(features, seqid, &qry_range, sty, err)))
      had_err = -1;
    if (!had_err && arguments->flattenfiles)
      gt_diagram_set_track_selector_func(d, flattened_file_track_selector,
                                         NULL);
    if (had_err || !(l = gt_layout_new(d, arguments->width, sty, err)))
      had_err = -1;
    if (!had_err)
      had_err = gt_layout_get_height(l, &height, err);
    if (!had_err) {
      ii = gt_image_info_new();

      if (strcmp(gt_str_get(arguments->format),"pdf")==0) {
        canvas = gt_canvas_cairo_file_new(sty, GT_GRAPHICS_PDF,
                                          arguments->width,
                                          height, ii, err);
      }
      else if (strcmp(gt_str_get(arguments->format),"ps")==0) {
        canvas = gt_canvas_cairo_file_new(sty, GT_GRAPHICS_PS,
                                          arguments->width,
                                          height, ii, err);
      }
      else if (strcmp(gt_str_get(arguments->format),"svg")==0) {
        canvas = gt_canvas_cairo_file_new(sty, GT_GRAPHICS_SVG,
                                          arguments->width,
                                          height, ii, err);
      }
      else {
        canvas = gt_canvas_cairo_file_new(sty, GT_GRAPHICS_PNG,
                                          arguments->width,
                                          height, ii, err);
      }
      if (!canvas)
        had_err = -1;
      if (!had_err) {
        had_err = gt_layout_sketch(l, canvas, err);
      }
      if (!had_err) {
        if (arguments->showrecmaps) {
          GtUword i;
          const GtRecMap *rm;
          for (i = 0; i < gt_image_info_num_of_rec_maps(ii) ;i++) {
            char buf[BUFSIZ];
            rm = gt_image_info_get_rec_map(ii, i);
            (void) gt_rec_map_format_html_imagemap_coords(rm, buf, BUFSIZ);
            printf("%s, %s\n",
                   buf,
                   gt_feature_node_get_type(gt_rec_map_get_genome_feature(rm)));
          }
        }
        if (arguments->use_streams) {
          GtFile *outfile;
          GtStr *str = gt_str_new();
          gt_canvas_cairo_file_to_stream((GtCanvasCairoFile*) canvas, str);
          outfile = gt_file_open(GT_FILE_MODE_UNCOMPRESSED, file, "w+", err);
          if (outfile) {
            gt_file_xwrite(outfile, gt_str_get_mem(str), gt_str_length(str));
            gt_file_delete(outfile);
          } else {
            had_err = -1;
          }
          gt_str_delete(str);
        } else {
          had_err = gt_canvas_cairo_file_to_file((GtCanvasCairoFile*) canvas,
                                                 file,
                                                 err);
        }
      }
    }
  }

  /* free */
  gt_free(seqid);
  gt_canvas_delete(canvas);
  gt_layout_delete(l);
  gt_image_info_delete(ii);
  gt_style_delete(sty);
  gt_diagram_delete(d);
  gt_array_delete(results);
  gt_str_delete(defaultstylefile);
  gt_feature_index_delete(features);

  return had_err;
}

GtTool* gt_sketch(void)
{
  return gt_tool_new(gt_sketch_arguments_new,
                     gt_sketch_arguments_delete,
                     gt_sketch_option_parser_new,
                     gt_sketch_arguments_check,
                     gt_sketch_runner);
}
