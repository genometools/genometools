/*
  Copyright (c) 2008-2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2009 Center for Bioinformatics, University of Hamburg

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

#include <time.h>
#include <string.h>
#include <cairo.h>
#if CAIRO_HAS_PS_SURFACE
#include <cairo-ps.h>
#endif
#if CAIRO_HAS_PDF_SURFACE
#include <cairo-pdf.h>
#endif
#if CAIRO_HAS_SVG_SURFACE
#include <cairo-svg.h>
#endif
#include "core/bioseq.h"
#include "core/cstr_api.h"
#include "core/fileutils_api.h"
#include "core/gtdatapath.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/option_api.h"
#include "core/str.h"
#include "core/unused_api.h"
#include "core/undef_api.h"
#include "core/versionfunc.h"
#include "core/warning_api.h"
#include "annotationsketch/canvas_cairo_context.h"
#include "annotationsketch/custom_track_gc_content_api.h"
#include "annotationsketch/diagram.h"
#include "extended/feature_index_memory.h"
#include "annotationsketch/gt_sketch_page.h"
#include "annotationsketch/layout.h"
#include "annotationsketch/style.h"
#include "annotationsketch/text_width_calculator.h"
#include "annotationsketch/text_width_calculator_cairo.h"

#define TEXT_SPACER         8  /* pt */
#define TIME_DATE_FORMAT    "%a, %b %d %Y - %T"

typedef struct {
  unsigned long width;
  double pwidth, pheight, theight;
  GtRange range;
  GtStr *seqid, *format, *stylefile, *text, *seqfile;
} SketchPageArguments;

static void *gt_sketch_page_arguments_new(void)
{
  SketchPageArguments *arguments = gt_malloc(sizeof *arguments);
  arguments->format = gt_str_new();
  arguments->seqid = gt_str_new();
  arguments->text = gt_str_new();
  arguments->stylefile = gt_str_new();
  arguments->seqfile = gt_str_new();
  arguments->range.start = arguments->range.end == GT_UNDEF_ULONG;
  return arguments;
}

static void gt_sketch_page_arguments_delete(void *tool_arguments)
{
  SketchPageArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->seqid);
  gt_str_delete(arguments->format);
  gt_str_delete(arguments->text);
  gt_str_delete(arguments->stylefile);
  gt_str_delete(arguments->seqfile);
  gt_free(arguments);
}

static GtOptionParser* gt_sketch_page_option_parser_new(void *tool_arguments)
{
  SketchPageArguments *arguments = tool_arguments;
  GtOptionParser *op;
  static const char *formats[] = {
#ifdef CAIRO_HAS_PDF_SURFACE
    "pdf",
#endif
#ifdef CAIRO_HAS_PS_SURFACE
    "ps",
#endif
    NULL
  };
  GtOption *o;
  op = gt_option_parser_new("outfile annotationfile",
                            "Draw a multi-page PDF/PS representation of "
                            "an annotation file.");
  o = gt_option_new_string("seqid", "sequence region to draw\n"
                                    "default: first in file",
                           arguments->seqid, NULL);
  gt_option_parser_add_option(op, o);
  gt_option_hide_default(o);

  o = gt_option_new_string("text", "text to show in header\n"
                                  "default: file name",
                           arguments->text, NULL);
  gt_option_parser_add_option(op, o);
  gt_option_hide_default(o);

  o = gt_option_new_double("fontsize", "header and footer font size "
                                       "(in points)",
                           &arguments->theight, 10.0);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_range("range", "range to draw (e.g. 100 10000)\n"
                                   "default: full range",
                          &arguments->range, NULL);
  gt_option_parser_add_option(op, o);
  gt_option_hide_default(o);

  o = gt_option_new_ulong_min("linewidth", "base width of a single "
                                           "repeated unit",
                              &arguments->width, 2000, 1000);
  gt_option_is_mandatory(o);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_double("width", "page width in millimeters "
                                    "(default: DIN A4)",
                           &arguments->pwidth, 210.0);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_double("height", "page height in millimeters "
                                     "(default: DIN A4)",
                           &arguments->pheight, 297.0);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_choice("format", "output format\n"
                                     "choose from: "
#ifdef CAIRO_HAS_PDF_SURFACE
                                       "pdf"
#ifdef CAIRO_HAS_PS_SURFACE
                                       "|"
#endif
#endif
#ifdef CAIRO_HAS_PS_SURFACE
                                       "ps"
#endif
                                       "",
                            arguments->format, formats[0], formats );
  gt_option_parser_add_option(op, o);

  o = gt_option_new_string("style", "style file to use\n"
                                    "default: gtdata/sketch/default.style",
                                arguments->stylefile,
                                gt_str_get(arguments->stylefile));
  gt_option_parser_add_option(op, o);
  gt_option_hide_default(o);

  o = gt_option_new_filename("seqfile", "sequence file for GC content view",
                                arguments->seqfile);
  gt_option_parser_add_option(op, o);
  gt_option_is_extended_option(o);

  gt_option_parser_set_min_max_args(op, 2, 2);
  return op;
}

static double mm_to_pt(double mm)
{
  return 2.8457598*mm;
}

static void draw_header(cairo_t *cr, const char *text, GT_UNUSED const char *fn,
                        const char *seqid, GT_UNUSED unsigned long pagenum,
                        GT_UNUSED double width, GT_UNUSED double height,
                        double theight)
{
  cairo_text_extents_t ext;
  time_t t;
  char buffer[BUFSIZ];
  struct tm *tmp;
  double xpos = TEXT_SPACER;
  cairo_save(cr);
  t = time(NULL);
  tmp = localtime(&t);
  cairo_set_font_size(cr, theight);
  if (tmp)
  {
    (void) strftime(buffer, BUFSIZ, TIME_DATE_FORMAT, tmp);
    cairo_text_extents(cr, buffer, &ext);
    cairo_move_to(cr, width - TEXT_SPACER - ext.width, TEXT_SPACER
                    + theight);
    cairo_set_source_rgba(cr, 0,0,0,1);
    cairo_show_text(cr, buffer);
  }
  cairo_text_extents(cr, text, &ext);
  cairo_move_to(cr, xpos, TEXT_SPACER + theight);
  cairo_set_source_rgba(cr, 0,0,0,1);
  cairo_show_text(cr, text);
  xpos += ext.width+3;
  cairo_move_to(cr, xpos, TEXT_SPACER + theight);
  cairo_set_source_rgba(cr, 0.7, 0.7, 0.7, 1);
  cairo_show_text(cr, ", sequence region: ");
  cairo_text_extents(cr, ", sequence region: ", &ext);
  xpos += ext.width + 10;
  cairo_set_source_rgba(cr, 0,0,0,1);
  cairo_show_text(cr, seqid);
  xpos = TEXT_SPACER;
  cairo_move_to(cr, xpos, height - 2*TEXT_SPACER - theight);
  (void) snprintf(buffer, BUFSIZ, "Page %lu", pagenum+1);
  cairo_show_text(cr, buffer);
  cairo_restore(cr);
}

static int gt_sketch_page_runner(GT_UNUSED int argc,
                                 const char **argv,
                                 int parsed_args,
                                 void *tool_arguments,
                                 GtError *err)
{
  SketchPageArguments *arguments = tool_arguments;
  int had_err = 0;
  GtFeatureIndex *features = NULL;
  GtRange qry_range, sequence_region_range;
  GtStyle *sty = NULL;
  GtStr *prog, *gt_style_file;
  GtDiagram *d = NULL;
  GtLayout *l = NULL;
  GtBioseq *bioseq = NULL;
  GtCanvas *canvas = NULL;
  char *seqid = NULL;
  const char *outfile = NULL;
  unsigned long start, height, num_pages = 0;
  double offsetpos, usable_height;
  cairo_surface_t *surf = NULL;
  cairo_t *cr = NULL;
  bool has_seqid;
  GtTextWidthCalculator *twc;
  gt_error_check(err);

  features = gt_feature_index_memory_new();

  if (cairo_version() < CAIRO_VERSION_ENCODE(1, 8, 6))
    gt_warning("Your cairo library (version %s) is older than version 1.8.6! "
               "These versions contain a bug which may result in "
               "corrupted PDF output!", cairo_version_string());

  /* get style */
  sty = gt_style_new(err);
  if (gt_str_length(arguments->stylefile) == 0)
  {
    prog = gt_str_new();
    gt_str_append_cstr_nt(prog, argv[0],
                          gt_cstr_length_up_to_char(argv[0], ' '));
    gt_style_file = gt_get_gtdata_path(gt_str_get(prog), err);
    gt_str_delete(prog);
    gt_str_append_cstr(gt_style_file, "/sketch/default.style");
  }
  else
  {
    gt_style_file = gt_str_ref(arguments->stylefile);
  }
  had_err = gt_style_load_file(sty, gt_str_get(gt_style_file), err);
  if (!had_err) {
    had_err = gt_feature_index_has_seqid(features, &has_seqid,
                                         gt_str_get(arguments->seqid), err);
  }

  outfile = argv[parsed_args];
  if (!had_err)
  {
    /* get features */
    had_err = gt_feature_index_add_gff3file(features, argv[parsed_args+1], err);
     if (!had_err && gt_str_length(arguments->seqid) == 0) {
      seqid = gt_feature_index_get_first_seqid(features, err);
      if (seqid == NULL)
      {
        gt_error_set(err, "GFF input file must contain a sequence region!");
        had_err = -1;
      }
    }
    else if (!had_err && !has_seqid)
    {
      gt_error_set(err, "sequence region '%s' does not exist in GFF input file",
                   gt_str_get(arguments->seqid));
      had_err = -1;
    }
    else if (!had_err)
      seqid = gt_str_get(arguments->seqid);
  }

  /* set text */
  if (gt_str_length(arguments->text) == 0)
  {
    gt_str_delete(arguments->text);
    arguments->text = gt_str_new_cstr(argv[parsed_args+1]);
  }

  if (!had_err)
  {
    /* set display range */
    had_err = gt_feature_index_get_range_for_seqid(features,
                                                   &sequence_region_range,
                                                   seqid, err);
  }
  if (!had_err)
  {
    qry_range.start = (arguments->range.start == GT_UNDEF_ULONG ?
                         sequence_region_range.start :
                         arguments->range.start);
    qry_range.end   = (arguments->range.end == GT_UNDEF_ULONG ?
                         sequence_region_range.end :
                         arguments->range.end);

    /* set output format */
    if (strcmp(gt_str_get(arguments->format), "pdf") == 0)
    {
      surf = cairo_pdf_surface_create(outfile,
                                      mm_to_pt(arguments->pwidth),
                                      mm_to_pt(arguments->pheight));
    }
    else if (strcmp(gt_str_get(arguments->format), "ps") == 0)
    {
      surf =  cairo_ps_surface_create(outfile,
                                      mm_to_pt(arguments->pwidth),
                                      mm_to_pt(arguments->pheight));
    }
    gt_log_log("created page with %.2f:%.2f dimensions\n",
                                                  mm_to_pt(arguments->pwidth),
                                                  mm_to_pt(arguments->pheight));

    offsetpos = TEXT_SPACER + arguments->theight + TEXT_SPACER;
    usable_height = mm_to_pt(arguments->pheight)
                              - arguments->theight
                              - arguments->theight
                              - 4*TEXT_SPACER;

    if (gt_str_length(arguments->seqfile) > 0) {
      bioseq = gt_bioseq_new(gt_str_get(arguments->seqfile), err);
    }

    cr = cairo_create(surf);
    cairo_set_font_size(cr, 8);
    twc = gt_text_width_calculator_cairo_new(cr, sty, err);
    for (start = qry_range.start; start <= qry_range.end;
         start += arguments->width)
    {
      GtRange single_range;
      GtCustomTrack *ct = NULL;
      const char *seq;
      single_range.start = start;
      single_range.end = start + arguments->width;

      if (had_err)
        break;

      d = gt_diagram_new(features, seqid, &single_range, sty, err);
      if (!d) {
        had_err = -1;
        break;
      }
      if (bioseq) {
        seq = gt_bioseq_get_sequence(bioseq, 0);
        ct = gt_custom_track_gc_content_new(seq,
                                      gt_bioseq_get_sequence_length(bioseq, 0),
                                      800, 70, 0.4, true);
        gt_diagram_add_custom_track(d, ct);
      }

      l = gt_layout_new_with_twc(d, mm_to_pt(arguments->width), sty, twc, err);
      had_err = gt_layout_get_height(l, &height, err);
      if (!had_err) {
        if (gt_double_smaller_double(usable_height - 10 - 2*TEXT_SPACER
              - arguments->theight, offsetpos + height))
        {
            draw_header(cr, gt_str_get(arguments->text), argv[parsed_args+1],
                        seqid, num_pages, mm_to_pt(arguments->pwidth),
                        mm_to_pt(arguments->pheight),
                        arguments->theight);
          cairo_show_page(cr);
          offsetpos = TEXT_SPACER + arguments->theight + TEXT_SPACER;
          num_pages++;
        }
        canvas = gt_canvas_cairo_context_new(sty,
                                             cr,
                                             offsetpos,
                                             mm_to_pt(arguments->pwidth),
                                             height,
                                             NULL,
                                             err);
        if (!canvas)
          had_err = -1;
        offsetpos += height;
        if (!had_err)
          had_err = gt_layout_sketch(l, canvas, err);
      }
      gt_canvas_delete(canvas);
      gt_layout_delete(l);
      gt_diagram_delete(d);
      if (ct)
        gt_custom_track_delete(ct);
    }
    draw_header(cr, gt_str_get(arguments->text), argv[parsed_args+1], seqid,
                num_pages, mm_to_pt(arguments->pwidth),
                mm_to_pt(arguments->pheight),
                arguments->theight);
    cairo_show_page(cr);
    num_pages++;
    gt_log_log("finished, should be %lu pages\n", num_pages);
    gt_text_width_calculator_delete(twc);
    cairo_destroy(cr);
    cairo_surface_flush(surf);
    cairo_surface_finish(surf);
    cairo_surface_destroy(surf);
    cairo_debug_reset_static_data();
    if (bioseq)
      gt_bioseq_delete(bioseq);
    gt_style_delete(sty);
    gt_free(seqid);
    gt_str_delete(gt_style_file);
    gt_feature_index_delete(features);
  }
  return had_err;
}

GtTool* gt_sketch_page(void)
{
  return gt_tool_new(gt_sketch_page_arguments_new,
                     gt_sketch_page_arguments_delete,
                     gt_sketch_page_option_parser_new,
                     NULL,
                     gt_sketch_page_runner);
}
