/*
  Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>

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
#include "core/cstr.h"
#include "core/fileutils.h"
#include "core/gtdatapath.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/option.h"
#include "core/unused_api.h"
#include "core/undef.h"
#include "core/versionfunc.h"
#include "annotationsketch/canvas_cairo_context.h"
#include "annotationsketch/custom_track_gc_content_api.h"
#include "annotationsketch/diagram.h"
#include "annotationsketch/feature_index_memory.h"
#include "annotationsketch/gt_sketch_page.h"
#include "annotationsketch/layout.h"
#include "annotationsketch/style.h"

#define HEADER_TEXT_HEIGHT 20  /* pt */
#define FOOTER_TEXT_HEIGHT 20  /* pt */
#define TEXT_SPACER         8  /* pt */
#define TIME_DATE_FORMAT    "%a, %b %d %Y - %T"

typedef struct {
  unsigned long width;
  double pwidth, pheight;
  GtRange range;
  GtStr *seqid, *format, *text;
} SketchPageArguments;

static void *gt_sketch_page_arguments_new(void)
{
  SketchPageArguments *arguments = gt_malloc(sizeof *arguments);
  arguments->format = gt_str_new();
  arguments->seqid = gt_str_new();
  arguments->text = gt_str_new();
  arguments->range.start = arguments->range.end == UNDEF_ULONG;
  return arguments;
}

static void gt_sketch_page_arguments_delete(void *tool_arguments)
{
  SketchPageArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->seqid);
  gt_str_delete(arguments->format);
  gt_str_delete(arguments->text);
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
  o = gt_option_new_string("seqid", "sequence region to draw",
                           arguments->seqid, NULL);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_string("text", "text to show in header",
                           arguments->text, NULL);
  gt_option_parser_add_option(op, o);

  o = gt_option_new_range("range", "range to draw",
                          &arguments->range, NULL);
  gt_option_parser_add_option(op, o);
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
  gt_option_parser_set_min_max_args(op, 2, 2);
  return op;
}

static inline double mm_to_pt(double mm)
{
  return 2.8457598*mm;
}

static void draw_header(cairo_t *cr, const char *text, GT_UNUSED const char *fn,
                        const char *seqid, GT_UNUSED unsigned long pagenum,
                        GT_UNUSED double width, GT_UNUSED double height)
{
  cairo_text_extents_t ext;
  time_t t;
  char buffer[BUFSIZ];
  struct tm *tmp;
  double xpos = TEXT_SPACER;
  cairo_save(cr);
  t = time(NULL);
  tmp = localtime(&t);
  cairo_set_font_size(cr, HEADER_TEXT_HEIGHT);
  if (tmp)
  {
    strftime(buffer, BUFSIZ, TIME_DATE_FORMAT, tmp);
    cairo_text_extents(cr, buffer, &ext);
    cairo_move_to(cr, width - TEXT_SPACER - ext.width, TEXT_SPACER + HEADER_TEXT_HEIGHT);
    cairo_set_source_rgba(cr, 0,0,0,1);
    cairo_show_text(cr, buffer);
  }
  cairo_text_extents(cr, text, &ext);
  cairo_move_to(cr, xpos, TEXT_SPACER + HEADER_TEXT_HEIGHT);
  cairo_set_source_rgba(cr, 0,0,0,1);
  cairo_show_text(cr, text);
  xpos += ext.width+3;
  cairo_move_to(cr, xpos, TEXT_SPACER + HEADER_TEXT_HEIGHT);
  cairo_set_source_rgba(cr, 0.7, 0.7, 0.7, 1);
  cairo_show_text(cr, ", sequence region: ");
  cairo_text_extents(cr, ", sequence region: ", &ext);
  xpos += ext.width + 10;
  cairo_move_to(cr, xpos, TEXT_SPACER + HEADER_TEXT_HEIGHT);
  cairo_set_source_rgba(cr, 0,0,0,1);
  cairo_show_text(cr, seqid);
  xpos = TEXT_SPACER;
  cairo_move_to(cr, xpos, height - 2*TEXT_SPACER - FOOTER_TEXT_HEIGHT);
  snprintf(buffer, BUFSIZ, "Page %lu", pagenum+1);
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
  const char *seqid = NULL, *outfile;
  unsigned long start, height, num_pages = 0;
  double offsetpos, usable_height;
  gt_error_check(err);
  cairo_surface_t *surf = NULL;
  cairo_t *cr = NULL;

  features = gt_feature_index_memory_new();

  /* get style */
  sty = gt_style_new(err);
  prog = gt_str_new();
  gt_str_append_cstr_nt(prog, argv[0],
                        gt_cstr_length_up_to_char(argv[0], ' '));
  gt_style_file = gt_get_gtdata_path(gt_str_get(prog), err);
  gt_str_delete(prog);
  gt_str_append_cstr(gt_style_file, "/sketch/default.style");
  gt_style_load_file(sty, gt_str_get(gt_style_file), err);

  /* get features */
  outfile = argv[parsed_args];
  had_err = gt_feature_index_add_gff3file(features, argv[parsed_args+1], err);
   if (!had_err && gt_str_length(arguments->seqid) == 0) {
    seqid = gt_feature_index_get_first_seqid(features);
    if (seqid == NULL)
    {
      gt_error_set(err, "GFF input file must contain a sequence region!");
      had_err = -1;
    }
  }
  else if (!had_err && !gt_feature_index_has_seqid(features,
                                                   gt_str_get(arguments->seqid)))
  {
    gt_error_set(err, "sequence region '%s' does not exist in GFF input file",
                 gt_str_get(arguments->seqid));
    had_err = -1;
  }
  else if (!had_err)
    seqid = gt_str_get(arguments->seqid);

  /* set text */
  if (gt_str_length(arguments->text) == 0)
  {
    gt_str_delete(arguments->text);
    arguments->text = gt_str_new_cstr(argv[parsed_args+1]);
  }

  /* set display range */
  gt_feature_index_get_range_for_seqid(features, &sequence_region_range, seqid);
  qry_range.start = (arguments->range.start == UNDEF_ULONG ?
                       sequence_region_range.start :
                       arguments->range.start);
  qry_range.end   = (arguments->range.end == UNDEF_ULONG ?
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

/*  bioseq = gt_bioseq_new("Drosophila_melanogaster"
 *                         ".BDGP5.4.51.dna.chromosome.2R.fa.gz", err); */
  offsetpos = TEXT_SPACER + HEADER_TEXT_HEIGHT + TEXT_SPACER;
  usable_height = mm_to_pt(arguments->pheight)
                            - HEADER_TEXT_HEIGHT
                            - FOOTER_TEXT_HEIGHT
                            - 4*TEXT_SPACER;

  /* do it! */
  cr = cairo_create(surf);
  cairo_set_font_size(cr, 8);
  for (start = qry_range.start; start <= qry_range.end;
       start += arguments->width)
  {
    GtRange single_range;
    /* GtCustomTrack *ct; */
    single_range.start = start;
    single_range.end = start + arguments->width;
  /* ct = gt_custom_track_gc_content_new(gt_bioseq_get_sequence(bioseq, 0),
                                        gt_bioseq_get_sequence_length(bioseq,
                                                                      0),
                                        300,
                                        20,
                                        0.365,
                                        true); */
    d = gt_diagram_new(features, seqid, &single_range, sty, err);
    gt_error_check(err);
    /* gt_diagram_add_custom_track(d, ct); */
    l = gt_layout_new(d, mm_to_pt(arguments->pwidth), sty, err);
    gt_error_check(err);
    height = gt_layout_get_height(l);
    if (offsetpos + height > usable_height - 10 - 2*TEXT_SPACER - FOOTER_TEXT_HEIGHT)
    {
        draw_header(cr, gt_str_get(arguments->text), argv[parsed_args+1], seqid, num_pages,
                    mm_to_pt(arguments->pwidth), mm_to_pt(arguments->pheight));
      cairo_show_page(cr);
      offsetpos = TEXT_SPACER + HEADER_TEXT_HEIGHT + TEXT_SPACER;
      num_pages++;
    }
    canvas = gt_canvas_cairo_context_new(sty,
                                         cr,
                                         offsetpos,
                                         mm_to_pt(arguments->pwidth),
                                         height,
                                         NULL);
    offsetpos += height;
    gt_layout_sketch(l, canvas, err);
    gt_error_check(err);
    /* gt_custom_track_delete(ct); */
    gt_canvas_delete(canvas);
    gt_layout_delete(l);
    gt_diagram_delete(d);
  }
  draw_header(cr, gt_str_get(arguments->text), argv[parsed_args+1], seqid, num_pages,
              mm_to_pt(arguments->pwidth), mm_to_pt(arguments->pheight));
  cairo_show_page(cr);
  num_pages++;
  gt_log_log("finished, should be %lu pages\n", num_pages);
  cairo_destroy(cr);
  cairo_surface_flush(surf);
  cairo_surface_finish(surf);
  cairo_surface_destroy(surf);
  cairo_debug_reset_static_data();
  if (bioseq)
    gt_bioseq_delete(bioseq);
  gt_style_delete(sty);
  gt_str_delete(gt_style_file);
  gt_feature_index_delete(features);
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
