/*
  Copyright (c) 2007-2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Gordon Gremme <gordon@gremme.org>,
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

#include <cairo.h>
#ifdef CAIRO_HAS_PDF_SURFACE
#include <cairo-pdf.h>
#endif
#ifdef CAIRO_HAS_PS_SURFACE
#include <cairo-ps.h>
#endif
#ifdef CAIRO_HAS_SVG_SURFACE
#include <cairo-svg.h>
#endif
#include <pango/pangocairo.h>

#include <math.h>
#include "annotationsketch/default_formats.h"
#include "annotationsketch/graphics_cairo_api.h"
#include "annotationsketch/graphics_rep.h"
#include "core/file.h"
#include "core/fileutils_api.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/str.h"
#include "core/unused_api.h"

struct GtGraphicsCairo {
  const GtGraphics parent_instance;
  cairo_t *cr;
  cairo_surface_t *surf;
  GtColor bg_color;
  GtStr *outbuf;
  GtGraphicsOutType type;
  double margin_x, margin_y, height, width;
  bool from_context;
  PangoLayout *layout;
  PangoFontDescription *desc;
  int font_height;
};

const GtGraphicsClass* gt_graphics_cairo_class(void);

#define gt_graphics_cairo_cast(G)\
        gt_graphics_cast(gt_graphics_cairo_class(), G)

static cairo_status_t str_write_func(void *closure, const unsigned char *data,
                                     unsigned int length)
{
  GtStr *stream = closure;
  gt_assert(stream);
  gt_str_append_cstr_nt(stream, (char*) data, length);
  return CAIRO_STATUS_SUCCESS;
}

/* to get crisp lines, round coordinates to .5 */
#define rnd_to_nhalf(num) (floor(num+0.5)+0.5)

void gt_graphics_cairo_set_font(GtGraphics *gg, const char *family,
                                FontSlant slant, FontWeight weight, double size)
{
  char buf[BUFSIZ];
  GtGraphicsCairo *g = gt_graphics_cairo_cast(gg);
  gt_assert(g && family && g->layout);

  snprintf(buf, BUFSIZ, "%s %s %s %d",
           family,
           (slant == SLANT_ITALIC) ? "Italic" : "",
           (weight == WEIGHT_BOLD) ? "Bold" : "",
           (int) size);
  g->desc = pango_font_description_from_string(buf);
  gt_assert(g->desc);
  pango_layout_set_font_description(g->layout, g->desc);
  pango_font_description_free(g->desc);
  g->font_height = (int) size;
}

void gt_graphics_cairo_initialize(GtGraphics *gg, GtGraphicsOutType type,
                                  unsigned int width, unsigned int height)
{
  GtGraphicsCairo *g = gt_graphics_cairo_cast(gg);
  g->outbuf = gt_str_new();
  switch (type)
  {
    case GT_GRAPHICS_PDF:
#ifdef CAIRO_HAS_PDF_SURFACE
      g->surf = cairo_pdf_surface_create_for_stream(str_write_func,
                                                    g->outbuf,
                                                    width,
                                                    height);
      break;
#endif
    case GT_GRAPHICS_PS:
#ifdef CAIRO_HAS_PS_SURFACE
      g->surf = cairo_ps_surface_create_for_stream(str_write_func,
                                                   g->outbuf,
                                                   width,
                                                   height);
      break;
#endif
    case GT_GRAPHICS_SVG:
#ifdef CAIRO_HAS_SVG_SURFACE
      g->surf = cairo_svg_surface_create_for_stream(str_write_func,
                                                    g->outbuf,
                                                    width,
                                                    height);
      break;
#endif
    case GT_GRAPHICS_PNG:
    default:
      g->surf = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
      break;
  }
  gt_assert(g->surf && cairo_surface_status(g->surf) == CAIRO_STATUS_SUCCESS);
  g->cr = cairo_create(g->surf);
  gt_assert(cairo_status(g->cr) == CAIRO_STATUS_SUCCESS);
  /* set background default to transparent */
  g->bg_color.red = g->bg_color.green  = 0.0;
  g->bg_color.blue = g->bg_color.alpha = 0.0;
  g->width = width;
  g->height = height;
  g->margin_x = g->margin_y = 20;
  cairo_set_line_join(g->cr, CAIRO_LINE_JOIN_ROUND);
  cairo_set_line_cap(g->cr, CAIRO_LINE_CAP_ROUND);
  g->type = type;
}

int gt_graphics_cairo_set_background_color(GtGraphics *gg, GtColor color)
{
  GtGraphicsCairo *g = gt_graphics_cairo_cast(gg);
  gt_assert(g);
  if (g->type != GT_GRAPHICS_PNG)
  {
    /* blending with a background is only supported with image output type */
    return -1;
  }
  else
  {
    g->bg_color.red = color.red;
    g->bg_color.green = color.green;
    g->bg_color.blue = color.blue;
    g->bg_color.alpha = color.alpha;
    return 0;
  }
}

void gt_graphics_cairo_draw_text(GtGraphics *gg, double x, double y,
                                 const char *text)
{
  PangoRectangle ink;
  GtGraphicsCairo *g = gt_graphics_cairo_cast(gg);
  gt_assert(g && text && g->layout);

  pango_layout_set_text(g->layout, text, -1);

  /* get text extents */
  pango_layout_get_pixel_extents(g->layout, &ink, NULL);

  if (gt_double_smaller_double(g->width, x+ink.width))
    return;

  cairo_set_source_rgb(g->cr, 0, 0, 0);
  cairo_move_to(g->cr, x, y-g->font_height);
  pango_cairo_show_layout(g->cr, g->layout);
}

void gt_graphics_cairo_draw_text_clip(GtGraphics *gg, double x, double y,
                                      const char *text)
{
  GtGraphicsCairo *g = gt_graphics_cairo_cast(gg);
  gt_assert(g && text && g->layout);

  pango_layout_set_text(g->layout, text, -1);

  cairo_save(g->cr);
  cairo_rectangle(g->cr,
                  g->margin_x,
                  g->margin_y,
                  g->width-2*g->margin_x,
                  g->height-2*g->margin_y);
  cairo_clip(g->cr);
  cairo_set_source_rgb(g->cr, 0, 0, 0);
  cairo_move_to(g->cr, x, y-g->font_height);
  pango_cairo_show_layout(g->cr, g->layout);
  cairo_restore(g->cr);
}

void gt_graphics_cairo_draw_text_centered(GtGraphics *gg, double x, double y,
                                          const char *text)
{
  PangoRectangle ink;
  GtGraphicsCairo *g = gt_graphics_cairo_cast(gg);
  gt_assert(g && text && g->layout);

  pango_layout_set_text(g->layout, text, -1);

  /* get text extents */
  pango_layout_get_pixel_extents(g->layout, &ink, NULL);

  cairo_set_source_rgb(g->cr, 0, 0, 0);
  /* draw text w/ its center at the given coords */
  cairo_move_to(g->cr, x-(ink.width/2)-1, y-g->font_height);
  pango_cairo_show_layout(g->cr, g->layout);
}

void gt_graphics_cairo_draw_text_right(GtGraphics *gg, double x, double y,
                                       const char *text)
{
  PangoRectangle ink;
  GtGraphicsCairo *g = gt_graphics_cairo_cast(gg);
  gt_assert(g && text && g->layout);

  pango_layout_set_text(g->layout, text, -1);

  /* get text extents */
  pango_layout_get_pixel_extents(g->layout, &ink, NULL);

  cairo_set_source_rgb(g->cr, 0, 0, 0);
  /* draw text w/ its right end at the given coords */
  cairo_move_to(g->cr, x-(ink.width)-1, y-g->font_height);
  pango_cairo_show_layout(g->cr, g->layout);
}

void gt_graphics_cairo_draw_colored_text(GtGraphics *gg, double x, double y,
                                         GtColor color, const char *text)
{
  GtGraphicsCairo *g = gt_graphics_cairo_cast(gg);
  gt_assert(g && text && g->layout);

  pango_layout_set_text(g->layout, text, -1);

  cairo_set_source_rgb(g->cr, color.red, color.green, color.blue);
  cairo_move_to(g->cr, x, y-g->font_height);
  pango_cairo_show_layout(g->cr, g->layout);
}

double gt_graphics_cairo_get_image_height(GtGraphics *gg)
{
  GtGraphicsCairo *g = gt_graphics_cairo_cast(gg);
  gt_assert(g);
  return g->height;
}

double gt_graphics_cairo_get_image_width(GtGraphics *gg)
{
  GtGraphicsCairo *g = gt_graphics_cairo_cast(gg);
  gt_assert(g);
  return g->width;
}

double gt_graphics_cairo_get_xmargins(GtGraphics *gg)
{
  GtGraphicsCairo *g = gt_graphics_cairo_cast(gg);
  gt_assert(g);
  return g->margin_x;
}

double gt_graphics_cairo_get_ymargins(GtGraphics *gg)
{
  GtGraphicsCairo *g = gt_graphics_cairo_cast(gg);
  gt_assert(g);
  return g->margin_y;
}

void gt_graphics_cairo_set_margins(GtGraphics *gg, double margin_x,
                                   double margin_y)
{
  GtGraphicsCairo *g = gt_graphics_cairo_cast(gg);
  gt_assert(g);
  g->margin_x = margin_x;
  g->margin_y = margin_y;
}

void gt_graphics_cairo_draw_line(GtGraphics *gg, double x, double y,
                                 double xto, double yto, GtColor color,
                                 double stroke_width)
{
  GtGraphicsCairo *g = gt_graphics_cairo_cast(gg);
  gt_assert(g);
  cairo_save(g->cr);
  cairo_move_to(g->cr, x, rnd_to_nhalf(y));
  cairo_line_to(g->cr, xto, rnd_to_nhalf(yto));
  cairo_set_line_width(g->cr, stroke_width);
  cairo_set_source_rgba(g->cr, color.red, color.green, color.blue, color.alpha);
  cairo_stroke(g->cr);
  cairo_restore(g->cr);
}

void gt_graphics_cairo_draw_horizontal_line(GtGraphics *gg, double x, double y,
                                            GtColor color, double width,
                                            double stroke_width)
{
  GtGraphicsCairo *g = gt_graphics_cairo_cast(gg);
  gt_assert(g);
  cairo_save(g->cr);
  cairo_set_line_width(g->cr, 1);
  cairo_move_to(g->cr, x, rnd_to_nhalf(y));
  cairo_rel_line_to(g->cr, width, 0);
  cairo_set_line_width(g->cr, stroke_width);
  cairo_set_source_rgba(g->cr, color.red, color.green, color.blue, color.alpha);
  cairo_stroke(g->cr);
  cairo_restore(g->cr);
}

void gt_graphics_cairo_draw_vertical_line(GtGraphics *gg, double x, double y,
                                          GtColor color, double length,
                                          double stroke_width)
{
  GtGraphicsCairo *g = gt_graphics_cairo_cast(gg);
  gt_assert(g);
  cairo_save(g->cr);
  cairo_move_to(g->cr, rnd_to_nhalf(x), rnd_to_nhalf(y));
  cairo_set_line_width(g->cr, 1);
  cairo_rel_line_to(g->cr, 0, floor(length));
  cairo_set_line_width(g->cr, stroke_width);
  cairo_set_source_rgba(g->cr, color.red, color.green, color.blue, color.alpha);
  cairo_stroke(g->cr);
  cairo_restore(g->cr);
}

double gt_graphics_cairo_get_text_width(GtGraphics *gg, const char* text)
{
  PangoRectangle rect;
  GtGraphicsCairo *g = gt_graphics_cairo_cast(gg);
  gt_assert(g && text && g->layout);

  pango_layout_set_text(g->layout, text, -1);
  /* get text extents */

  pango_layout_get_pixel_extents(g->layout, &rect, NULL);

  gt_assert(gt_double_smaller_double(0, rect.width));
  return rect.width;
}

double gt_graphics_cairo_get_text_height(GtGraphics *gg)
{
  GtGraphicsCairo *g = gt_graphics_cairo_cast(gg);
  gt_assert(g);

  return g->font_height;
}

void gt_graphics_cairo_draw_box(GtGraphics *gg, double x, double y,
                                double width, double height,
                                GtColor fill_color, ArrowStatus arrow_status,
                                double arrow_width, double stroke_width,
                                GtColor stroke_color, bool dashed)
{
  GtGraphicsCairo *g = gt_graphics_cairo_cast(gg);
  double dashes[]={2.0};
  bool widthdiff_geq0;
  gt_assert(g);

  /* save cairo context */
  cairo_save(g->cr);
  cairo_rectangle(g->cr, rnd_to_nhalf(g->margin_x), g->margin_y,
                  rnd_to_nhalf(g->width-2*g->margin_x),
                  g->height-2*g->margin_y);
  cairo_clip(g->cr);

  widthdiff_geq0 = gt_double_smaller_double(0, width - arrow_width);
  /* construct shape of the box or arrow */
  switch (arrow_status)
  {
    case ARROW_RIGHT:
      cairo_move_to(g->cr, rnd_to_nhalf(x), rnd_to_nhalf(y));
      if (widthdiff_geq0)
        cairo_line_to(g->cr, x + width - arrow_width, rnd_to_nhalf(y));
      cairo_line_to(g->cr, rnd_to_nhalf(x + width),
                           rnd_to_nhalf(y + height / 2));
      if (widthdiff_geq0)
        cairo_line_to(g->cr, x + width - arrow_width, rnd_to_nhalf(y + height));
      cairo_line_to(g->cr, rnd_to_nhalf(x), rnd_to_nhalf(y + height));
      cairo_close_path(g->cr);
      break;
    case ARROW_LEFT:
      cairo_move_to(g->cr, rnd_to_nhalf(x + width), rnd_to_nhalf(y));
      if (widthdiff_geq0) {
        cairo_line_to(g->cr, rnd_to_nhalf(x + arrow_width), rnd_to_nhalf(y));
      }
      cairo_line_to(g->cr, rnd_to_nhalf(x), rnd_to_nhalf(y + height / 2));
      if (widthdiff_geq0) {
        cairo_line_to(g->cr, rnd_to_nhalf(x + arrow_width),
                             rnd_to_nhalf(y + height));
      }
      cairo_line_to(g->cr, rnd_to_nhalf(x + width), rnd_to_nhalf(y + height));
      cairo_close_path(g->cr);
      break;
    case ARROW_BOTH:
      cairo_move_to(g->cr, rnd_to_nhalf(x), rnd_to_nhalf(y + height/2));
      if (gt_double_smaller_double(width, 2*arrow_width))
      {
        cairo_line_to(g->cr, rnd_to_nhalf(x + width/2), rnd_to_nhalf(y));
        cairo_line_to(g->cr, rnd_to_nhalf(x + width),
                             rnd_to_nhalf(y + height/2));
        cairo_line_to(g->cr, rnd_to_nhalf(x + width/2),
                             rnd_to_nhalf(y + height));
      }
      else
      {
        cairo_line_to(g->cr, rnd_to_nhalf(x + arrow_width), rnd_to_nhalf(y));
        cairo_line_to(g->cr, rnd_to_nhalf(x + width - arrow_width),
                             rnd_to_nhalf(y));
        cairo_line_to(g->cr, rnd_to_nhalf(x + width),
                             rnd_to_nhalf(y + height/2));
        cairo_line_to(g->cr, rnd_to_nhalf(x + width - arrow_width),
                             rnd_to_nhalf(y + height));
        cairo_line_to(g->cr, rnd_to_nhalf(x + arrow_width), y + height);
      }
      cairo_close_path(g->cr);
      break;
    case ARROW_NONE:
      cairo_rectangle(g->cr, rnd_to_nhalf(x), rnd_to_nhalf(y), width, height);
   }
   /* fill area */
   cairo_set_source_rgba(g->cr, fill_color.red,
                                fill_color.green,
                                fill_color.blue,
                                fill_color.alpha);
   cairo_fill_preserve(g->cr);
   /* draw outline */
   cairo_set_line_width(g->cr, stroke_width);
   cairo_set_source_rgba(g->cr, stroke_color.red,
                                stroke_color.green,
                                stroke_color.blue,
                                stroke_color.alpha);
   if (dashed)
     cairo_set_dash(g->cr, dashes, 1, (double) 0);
   cairo_stroke(g->cr);
   /* restore cairo context */
   cairo_restore(g->cr);
}

void gt_graphics_cairo_draw_dashes(GtGraphics *gg, double x, double y,
                                   double width, double height,
                                   ArrowStatus arrow_status, double arrow_width,
                                   double stroke_width, GtColor stroke_color)
{
  GtGraphicsCairo *g = gt_graphics_cairo_cast(gg);
  double dashes[] = {2.0};
  gt_assert(g);
  /* save cairo context */
  cairo_save(g->cr);
  cairo_rectangle(g->cr, g->margin_x, g->margin_y,
                  g->width-2*g->margin_x, g->height-2*g->margin_y);
  cairo_clip(g->cr);
  cairo_set_line_width(g->cr, stroke_width);
  cairo_set_source_rgba(g->cr, stroke_color.red,
                               stroke_color.green,
                               stroke_color.blue,
                               stroke_color.alpha);
  /* create arrowhead path */
  if (gt_double_smaller_double(0, width - arrow_width))
  {
    switch (arrow_status)
    {
      case ARROW_RIGHT:
        cairo_move_to(g->cr, rnd_to_nhalf(x + width - arrow_width), y);
        cairo_line_to(g->cr, rnd_to_nhalf(x + width), y + height / 2);
        cairo_line_to(g->cr, rnd_to_nhalf(x + width - arrow_width), y + height);
        /* draw arrowhead */
        cairo_stroke(g->cr);
        break;
      case ARROW_LEFT:
        cairo_move_to(g->cr, rnd_to_nhalf(x + arrow_width), y);
        cairo_line_to(g->cr, rnd_to_nhalf(x), y + height / 2);
        cairo_line_to(g->cr, rnd_to_nhalf(x + arrow_width), y + height);
        /* draw arrowhead */
        cairo_stroke(g->cr);
        break;
      case ARROW_BOTH:
        cairo_move_to(g->cr, rnd_to_nhalf(x + width - arrow_width), y);
        cairo_line_to(g->cr, rnd_to_nhalf(x + width), y + height / 2);
        cairo_line_to(g->cr, rnd_to_nhalf(x + width - arrow_width), y + height);
        cairo_move_to(g->cr, rnd_to_nhalf(x + arrow_width), y);
        cairo_line_to(g->cr, rnd_to_nhalf(x), y + height / 2);
        cairo_line_to(g->cr, rnd_to_nhalf(x + arrow_width), y + height);
        cairo_stroke(g->cr);
        break;
      case ARROW_NONE: break;
    }
  }
  /* draw dashes */
  cairo_set_dash(g->cr, dashes, 1, (double) 0);
  cairo_move_to(g->cr, rnd_to_nhalf(x), rnd_to_nhalf(y+(height/2)));
  cairo_rel_line_to(g->cr, rnd_to_nhalf(width), 0);
  cairo_stroke(g->cr);
  /* restore cairo context */
  cairo_restore(g->cr);
}

void gt_graphics_cairo_draw_caret(GtGraphics *gg, double x, double y,
                                  double width,double height,
                                  ArrowStatus arrow_status, double arrow_width,
                                  double stroke_width, GtColor stroke_color)
{
  GtGraphicsCairo *g = gt_graphics_cairo_cast(gg);
  gt_assert(g);

  /* save cairo context */
  cairo_save(g->cr);
  cairo_rectangle(g->cr, g->margin_x, g->margin_y,
                  g->width-2*g->margin_x, g->height-2*g->margin_y);
  cairo_clip(g->cr);
  /* create caret path */
  switch (arrow_status)
  {
    case ARROW_RIGHT:
      if (gt_double_smaller_double(0, width - arrow_width))
      {
        width -= arrow_width;
        cairo_move_to(g->cr, x, y + (height/2));
        cairo_line_to(g->cr, x + (width/2), y);
        cairo_line_to(g->cr, x + width, y + (height/2));
        cairo_line_to(g->cr, x + width + arrow_width, y + (height/2));
        /* add arrowhead */
        cairo_move_to(g->cr, x + width, y);
        cairo_line_to(g->cr, x + width + arrow_width, y + height / 2);
        cairo_line_to(g->cr, x + width, y + height);
      }
      break;
    case ARROW_LEFT:
      if (gt_double_smaller_double(0, width - arrow_width))
      {
        width -= arrow_width;
        /* draw arrowhead */
        cairo_move_to(g->cr, x + arrow_width , y);
        cairo_line_to(g->cr, x, y + (height/2));
        cairo_line_to(g->cr, x + arrow_width, y + height);
        cairo_move_to(g->cr, x, y + (height/2));
        cairo_line_to(g->cr, x + arrow_width, y + (height/2));
        cairo_line_to(g->cr, x + arrow_width + (width/2), y);
        cairo_line_to(g->cr, x + arrow_width + width, y + (height/2));
      }
      break;
    case ARROW_BOTH: /* XXX */
    case ARROW_NONE:
     cairo_move_to(g->cr, x, y + (height/2));
     cairo_line_to(g->cr, x + (width/2), y);
     cairo_line_to(g->cr, x + width, y + (height/2));
   }
   /* draw caret */
   cairo_set_line_width(g->cr, stroke_width);
   cairo_set_source_rgba(g->cr, stroke_color.red,
                                stroke_color.green,
                                stroke_color.blue,
                                stroke_color.alpha);
   cairo_stroke(g->cr);
   /* restore cairo context */
   cairo_restore(g->cr);
}

void gt_graphics_cairo_draw_rectangle(GtGraphics *gg, double x, double y,
                                      bool filled, GtColor fill_color,
                                      bool outlined, GtColor outline_color,
                                      double outline_width, double width,
                                      double height)
{
  GtGraphicsCairo *g = gt_graphics_cairo_cast(gg);
  gt_assert(g);

  /* save cairo context */
  cairo_save(g->cr);
  cairo_new_path(g->cr);
  cairo_rectangle(g->cr, x, y, width, height);
  if (filled) {
    cairo_set_source_rgba(g->cr, fill_color.red,
                                 fill_color.green,
                                 fill_color.blue,
                                 fill_color.alpha);
    cairo_fill_preserve(g->cr);
  }
  if (outlined) {
    cairo_set_line_width(g->cr, outline_width);
    cairo_set_source_rgba(g->cr, outline_color.red,
                                 outline_color.green,
                                 outline_color.blue,
                                 outline_color.alpha);
    cairo_stroke(g->cr);
  }
  /* restore cairo context */
  cairo_restore(g->cr);
}

void gt_graphics_cairo_draw_arrowhead(GtGraphics *gg, double x, double y,
                                      GtColor color, ArrowStatus arrow_status)
{
  GtGraphicsCairo *g = gt_graphics_cairo_cast(gg);
  double arrow_height = 8, arrow_width = 5;
  gt_assert(g);
  /* save cairo context */
  cairo_save(g->cr);
  cairo_reset_clip(g->cr);
  cairo_set_source_rgb(g->cr, color.red, color.green, color.blue);
  switch (arrow_status)
  {
    case ARROW_LEFT:
      cairo_move_to(g->cr, x+arrow_width, y);
      cairo_line_to(g->cr, x, y+(arrow_height/2));
      cairo_line_to(g->cr, x+arrow_width, y+arrow_height);
      cairo_close_path(g->cr);
      /* fill area */
      cairo_fill_preserve(g->cr);
      cairo_stroke(g->cr);
      break;
    case ARROW_RIGHT:
      cairo_move_to(g->cr, x, y);
      cairo_line_to(g->cr, x+arrow_width, y+(arrow_height/2));
      cairo_line_to(g->cr, x, y+arrow_height);
      cairo_close_path(g->cr);
      /* fill area */
      cairo_fill_preserve(g->cr);
      cairo_stroke(g->cr);
      break;
    case ARROW_BOTH: /* XXX */
    case ARROW_NONE: break;
  }
  /* restore cairo context */
  cairo_restore(g->cr);
}

int gt_graphics_cairo_save_to_file(const GtGraphics *gg, const char *filename,
                                GtError *err)
{
  const GtGraphicsCairo *g = (const GtGraphicsCairo*) gg;
  cairo_surface_t *bgsurf = NULL;
  cairo_t *bgc = NULL;
  cairo_status_t rval;
  GtFile *outfile;
  gt_error_check(err);
  gt_assert(g && filename);

  /* do nothing if no surface was created */
  if (g->from_context)
    return 0;
  switch (g->type)
  {
    case GT_GRAPHICS_PNG:
      /* blend rendered image with background color */
      bgsurf = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, g->width,
                                          g->height);
      bgc = cairo_create(bgsurf);
      cairo_set_source_rgba(bgc, g->bg_color.red, g->bg_color.green,
                                 g->bg_color.blue, g->bg_color.alpha);
      cairo_paint(bgc);
      cairo_set_source_surface(bgc, g->surf, 0, 0);
      cairo_paint(bgc);
      rval = cairo_surface_write_to_png(bgsurf, filename);
      gt_assert(rval == CAIRO_STATUS_SUCCESS ||
                rval == CAIRO_STATUS_WRITE_ERROR);
      if (rval == CAIRO_STATUS_WRITE_ERROR)
      {
        cairo_destroy(bgc);
        cairo_surface_destroy(bgsurf);
        gt_error_set(err, "an I/O error occurred while attempting "
                          "to write image file \"%s\"", filename);
        return -1;
      }
      cairo_destroy(bgc);
      cairo_surface_destroy(bgsurf);
      break;
    default:
      cairo_show_page(g->cr);
      cairo_surface_flush(g->surf);
      cairo_surface_finish(g->surf);
      outfile = gt_file_open(GT_FILE_MODE_UNCOMPRESSED, filename, "w+", err);
      if (outfile)
      {
        gt_file_xwrite(outfile, gt_str_get_mem(g->outbuf),
                       gt_str_length(g->outbuf));
        gt_file_delete(outfile);
      } else return -1;
      break;
  }
  return 0;
}

void gt_graphics_cairo_save_to_stream(const GtGraphics *gg, GtStr *stream)
{
  const GtGraphicsCairo *g = (const GtGraphicsCairo*) gg;
  GT_UNUSED cairo_status_t rval;
  cairo_surface_t *bgsurf = NULL;
  cairo_t *bgc = NULL;
  gt_assert(g && stream);

  /* do nothing if no surface was created */
  if (g->from_context)
    return;
  switch (g->type)
  {
    case GT_GRAPHICS_PNG:
      /* blend rendered image with background color */
      bgsurf = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, g->width,
                                          g->height);
      bgc = cairo_create(bgsurf);
      cairo_set_source_rgba(bgc, g->bg_color.red, g->bg_color.green,
                                 g->bg_color.blue, g->bg_color.alpha);
      cairo_paint(bgc);
      cairo_set_source_surface(bgc, g->surf, 0, 0);
      cairo_paint(bgc);
      rval = cairo_surface_write_to_png_stream(bgsurf, str_write_func, stream);
      gt_assert(rval == CAIRO_STATUS_SUCCESS); /* str_write_func() is sane */
      cairo_destroy(bgc);
      cairo_surface_destroy(bgsurf);
      break;
    default:
      cairo_show_page(g->cr);
      cairo_surface_flush(g->surf);
      cairo_surface_finish(g->surf);
      gt_str_reset(stream);
      gt_str_append_str(stream, g->outbuf);
      gt_assert(gt_str_length(stream) > 0);
  }
}

void gt_graphics_cairo_delete(GtGraphics *gg)
{
  GtGraphicsCairo *g;
  if (!gg) return;
  g = (GtGraphicsCairo*) gg;
  if (!g->from_context) /* do not attempt to destroy foreign contexts */
    cairo_destroy(g->cr);
  if (g->surf)
    cairo_surface_destroy(g->surf); /* reference counted */
  if (g->outbuf)
    gt_str_delete(g->outbuf);
  g_object_unref(g->layout);
}

const GtGraphicsClass* gt_graphics_cairo_class(void)
{
  static const GtGraphicsClass *gc = NULL;
  if (!gc)
  {
    gc = gt_graphics_class_new(sizeof (GtGraphicsCairo),
                               gt_graphics_cairo_draw_text,
                               gt_graphics_cairo_draw_text_clip,
                               gt_graphics_cairo_draw_text_centered,
                               gt_graphics_cairo_draw_text_right,
                               gt_graphics_cairo_draw_colored_text,
                               gt_graphics_cairo_get_text_height,
                               gt_graphics_cairo_get_text_width,
                               gt_graphics_cairo_set_background_color,
                               gt_graphics_cairo_set_font,
                               gt_graphics_cairo_get_image_width,
                               gt_graphics_cairo_get_image_height,
                               gt_graphics_cairo_get_xmargins,
                               gt_graphics_cairo_get_ymargins,
                               gt_graphics_cairo_set_margins,
                               gt_graphics_cairo_draw_line,
                               gt_graphics_cairo_draw_horizontal_line,
                               gt_graphics_cairo_draw_vertical_line,
                               gt_graphics_cairo_draw_box,
                               gt_graphics_cairo_draw_dashes,
                               gt_graphics_cairo_draw_caret,
                               gt_graphics_cairo_draw_rectangle,
                               gt_graphics_cairo_draw_arrowhead,
                               gt_graphics_cairo_draw_curve_data,
                               gt_graphics_cairo_save_to_file,
                               gt_graphics_cairo_save_to_stream,
                               gt_graphics_cairo_delete);
  }
  return gc;
}

GtGraphics* gt_graphics_cairo_new(GtGraphicsOutType type,
                                  unsigned int width,
                                  unsigned int height)
{
  GtGraphics *g;
  GtGraphicsCairo *gc;
  char buf[64];
  g = gt_graphics_create(gt_graphics_cairo_class());
  gc = gt_graphics_cairo_cast(g);
  gt_graphics_cairo_initialize(g, type, width, height);
  gc->layout = pango_cairo_create_layout(gc->cr);
  pango_layout_set_width(gc->layout, -1);
  gt_assert(gc->layout);
  snprintf(buf, 64, "Sans %d", TEXT_SIZE_DEFAULT);
  gc->desc = pango_font_description_from_string(buf);
  pango_layout_set_font_description(gc->layout, gc->desc);
  pango_font_description_free(gc->desc);
  return g;
}

GtGraphics* gt_graphics_cairo_new_from_context(cairo_t *context,
                                               unsigned int width,
                                               unsigned int height)
{
  GtGraphics *g;
  GtGraphicsCairo *gc;
  char buf[64];
  g = gt_graphics_create(gt_graphics_cairo_class());
  gc = gt_graphics_cairo_cast(g);
  gc->width = width;
  gc->height = height;
  gc->margin_x = gc->margin_y = 20;
  gc->from_context = true;
  gc->cr = context;
  gc->layout = pango_cairo_create_layout(gc->cr);
  pango_layout_set_width(gc->layout, -1);
  gt_assert(gc->layout);
  snprintf(buf, 64, "Sans %d", TEXT_SIZE_DEFAULT);
  gc->desc = pango_font_description_from_string(buf);
  pango_layout_set_font_description(gc->layout, gc->desc);
  pango_font_description_free(gc->desc);
  cairo_set_line_join(context, CAIRO_LINE_JOIN_ROUND);
  cairo_set_line_cap(context, CAIRO_LINE_CAP_ROUND);
  return g;
}

void gt_graphics_cairo_draw_curve_data(GtGraphics *gg, double x, double y,
                                       GtColor color, double data[],
                                       GtUword ndata, GtRange valrange,
                                       GtUword height)
{
  GtUword i, rnglen;
  double xpos;
  GtGraphicsCairo *g;
  g = gt_graphics_cairo_cast(gg);

  xpos = (((double)g->width-2*g->margin_x)/((double)ndata-1));
  rnglen = valrange.end - valrange.start;
  cairo_save(g->cr);
  cairo_move_to(g->cr,
                x,
                y+(1-(data[0]-valrange.start)/(double) rnglen)*height);
  for (i=1;i<ndata;i++)
  {
    double val, pval;
    if (gt_double_smaller_double(data[i], valrange.start)
          || gt_double_smaller_double(valrange.end, data[i]))
      break;
    val  = (double) (data[i]-valrange.start)/(double) rnglen;
    pval = (double) (data[i-1]-valrange.start)/(double) rnglen;
    gt_assert(val <= 1 && val >= 0 && pval >= 0 && pval <= 1);
    cairo_curve_to(g->cr,
                   x + ((i-0.5) * xpos),
                   y + (1-pval)   * height,
                   x + ((i-0.5) * xpos),
                   y + (1-val)    * height,
                   x + (i       * xpos),
                   y + (1-val)    * height);
  }
  cairo_set_source_rgba(g->cr, color.red,
                               color.green,
                               color.blue,
                               color.alpha);
  cairo_stroke(g->cr);
  cairo_restore(g->cr);
}
