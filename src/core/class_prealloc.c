/*
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#include "core/class_prealloc.h"
#include "extended/comment_node_api.h"
#include "extended/feature_node.h"
#include "extended/region_node_api.h"
#include "extended/sequence_node_api.h"
#include "extended/visitor_stream.h"
#ifndef WITHOUT_CAIRO
#include "annotationsketch/feature_index_memory.h"
#include "annotationsketch/feature_visitor.h"
#include "annotationsketch/canvas_cairo_file.h"
#include "annotationsketch/canvas_cairo_context.h"
#include "annotationsketch/custom_track_example.h"
#include "annotationsketch/custom_track_gc_content.h"
#include "annotationsketch/custom_track_script_wrapper.h"
#include "annotationsketch/text_width_calculator_cairo.h"
#include "annotationsketch/line_breaker_captions.h"
#include "annotationsketch/line_breaker_bases.h"
#include "annotationsketch/graphics_cairo.h"
#endif

void gt_class_prealloc_run(void)
{
  (void) gt_feature_node_class();
  (void) gt_comment_node_class();
  (void) gt_region_node_class();
  (void) gt_sequence_node_class();
  (void) gt_visitor_stream_class();
#ifndef WITHOUT_CAIRO
  (void) gt_feature_visitor_class();
  (void) gt_feature_index_memory_class();
  (void) gt_canvas_cairo_context_class();
  (void) gt_canvas_cairo_file_class();
  (void) gt_custom_track_example_class();
  (void) gt_custom_track_gc_content_class();
  (void) gt_custom_track_script_wrapper_class();
  (void) gt_text_width_calculator_cairo_class();
  (void) gt_line_breaker_bases_class();
  (void) gt_line_breaker_captions_class();
  (void) gt_graphics_cairo_class();
#endif
}
