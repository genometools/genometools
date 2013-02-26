/*
  Copyright (c) 2008-2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2012 Center for Bioinformatics, University of Hamburg

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
#include <pango/pangocairo.h>
#include "core/class_alloc_lock.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/unused_api.h"
#include "annotationsketch/default_formats.h"
#include "annotationsketch/style.h"
#include "annotationsketch/text_width_calculator.h"
#include "annotationsketch/text_width_calculator_cairo.h"
#include "annotationsketch/text_width_calculator_rep.h"

#define GT_TWC_FORMAT CAIRO_FORMAT_ARGB32
#define GT_TWC_WIDTH  500
#define GT_TWC_HEIGHT 60

struct GtTextWidthCalculatorCairo {
  const GtTextWidthCalculator parent_instance;
  GtStyle *style;
  cairo_t *context;
  cairo_surface_t *mysurf;
  PangoLayout *layout;
  PangoFontDescription *desc;
  bool own_context;
};

#define gt_text_width_calculator_cairo_cast(TWC)\
        gt_text_width_calculator_cast(gt_text_width_calculator_cairo_class(),\
                                      TWC)

double gt_text_width_calculator_cairo_get_text_width(GtTextWidthCalculator *twc,
                                                     const char *text,
                                                     GT_UNUSED GtError *err)
{
  GtTextWidthCalculatorCairo *twcc;
  PangoRectangle rect;
  gt_assert(twc && text);
  twcc = gt_text_width_calculator_cairo_cast(twc);

  /* redo layout */
  pango_layout_set_text(twcc->layout, text, -1);

  /* get extents */
  pango_layout_get_pixel_extents(twcc->layout, &rect, NULL);

  if (twcc->style)
    cairo_restore(twcc->context);
  gt_assert(gt_double_smaller_double(0, rect.width));
  return rect.width;
}

void gt_text_width_calculator_cairo_delete(GtTextWidthCalculator *twc)
{
  GtTextWidthCalculatorCairo *twcc;
  if (!twc) return;
  twcc = gt_text_width_calculator_cairo_cast(twc);
  g_object_unref(twcc->layout);
  if (twcc->style)
    gt_style_delete(twcc->style);
  if (twcc->own_context)
  {
    cairo_destroy(twcc->context);
    cairo_surface_destroy(twcc->mysurf);
  }
}

const GtTextWidthCalculatorClass* gt_text_width_calculator_cairo_class(void)
{
  static const GtTextWidthCalculatorClass *twcc = NULL;
  gt_class_alloc_lock_enter();
  if (!twcc)
  {
    twcc = gt_text_width_calculator_class_new(
                                  sizeof (GtTextWidthCalculatorCairo),
                                  gt_text_width_calculator_cairo_get_text_width,
                                  gt_text_width_calculator_cairo_delete);
  }
  gt_class_alloc_lock_leave();
  return twcc;
}

GtTextWidthCalculator* gt_text_width_calculator_cairo_new(cairo_t *context,
                                                          GtStyle *style,
                                                          GtError *err)
{
  GtTextWidthCalculatorCairo *twcc;
  GtTextWidthCalculator *twc;
  double theight = TEXT_SIZE_DEFAULT;
  char buf[64];
  twc = gt_text_width_calculator_create(gt_text_width_calculator_cairo_class());
  twcc = gt_text_width_calculator_cairo_cast(twc);
  if (style)
    twcc->style = gt_style_ref(style);
  if (!context)
  {
    twcc->mysurf = cairo_image_surface_create(GT_TWC_FORMAT, GT_TWC_WIDTH,
                                              GT_TWC_HEIGHT);
    twcc->context = cairo_create(twcc->mysurf);
    twcc->own_context = true;
  } else {
    twcc->context = context;
    twcc->own_context = false;
  }
  if (twcc->style)
  {
    if (gt_style_get_num(twcc->style,
                         "format", "block_caption_font_size",
                         &theight, NULL, err) == GT_STYLE_QUERY_ERROR) {
      gt_text_width_calculator_delete(twc);
      return NULL;
    }
    cairo_save(twcc->context);
  }
  twcc->layout = pango_cairo_create_layout(twcc->context);
  snprintf(buf, 64, "Sans %d", (int) theight);
  twcc->desc = pango_font_description_from_string(buf);
  pango_layout_set_font_description(twcc->layout, twcc->desc);
  pango_font_description_free(twcc->desc);
  return twc;
}
