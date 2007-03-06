/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg

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

#include <assert.h>
#include <cairo.h>
#include "cairo_visitor.h"
#include "genome_visitor_rep.h"
#include "minmax.h"
#include "undef.h"
#include "xansi.h"

#define TRACK_HEIGHT 		50
#define SPACE        		.05
#define ROOM         		.9
#define TEXT_POSITION   	20
#define FEATURE_POSITION	30
#define EXON_HEIGHT		10
#define EXON_ARROW_WIDTH        8

struct _Cairo_visitor {
  const Genome_visitor parent_instance;
  int width,
      height,
      global_track_number;
  unsigned int number_of_tracks;
  unsigned long from,
                to;
  cairo_t *cr;
  cairo_surface_t *surf;
  char *png_filename;
  Range drawed_sequence_range,
        last_range;
  unsigned int drawed_sequence_range_is_defined : 1,
               last_range_is_defined            : 1;
};

typedef struct {
  Cairo_visitor *cv;
  unsigned int children_overlap : 1;
  int local_track_number;
} Show_children_info;

#define cairo_visitor_cast(GV)\
        genome_visitor_cast(cairo_visitor_class(), GV)

static void cairo_visitor_free(Genome_visitor *gv)
{
  Cairo_visitor *cv = cairo_visitor_cast(gv);
  assert(cv->png_filename);
  assert(cv->width); /* the width has to be positive */
  assert(cv->height); /* the height has to be positive */
  (void) cairo_surface_write_to_png(cv->surf, cv->png_filename);
  cairo_surface_destroy(cv->surf); /* reference counted */
  /* we check this after writing the png to simplify debugging */
  assert(cv->global_track_number <= cv->number_of_tracks);
  cairo_destroy(cv->cr);
}

static void draw_exon_box(Cairo_visitor *cv, Genome_feature *gf, double width,
                          int track_number)
{
  Range feature_range = genome_node_get_range((Genome_node*) gf);
  Strand feature_strand = genome_feature_get_strand(gf);
  double x, y, height;

  x =  cv->width * SPACE + cv->width * ROOM *
       ((double) (feature_range.start - cv->drawed_sequence_range.start + 1) /
        (double) range_length(cv->drawed_sequence_range));
  y = track_number * TRACK_HEIGHT + FEATURE_POSITION;
  height = EXON_HEIGHT;

  cairo_set_source_rgb(cv->cr, 0, 0, 1);
  switch (feature_strand) {
    case STRAND_FORWARD:
      cairo_move_to(cv->cr, x, y);
      if (width - EXON_ARROW_WIDTH > 0)
        cairo_rel_line_to(cv->cr, width - EXON_ARROW_WIDTH, 0);
      cairo_line_to(cv->cr, x + width, y + height / 2);
      if (width - EXON_ARROW_WIDTH > 0)
        cairo_line_to(cv->cr, x + width - EXON_ARROW_WIDTH, y + height);
      cairo_line_to(cv->cr, x, y + height);
      cairo_close_path(cv->cr);
      break;
    case STRAND_REVERSE:
      cairo_move_to(cv->cr, x + width, y);
      if (width - EXON_ARROW_WIDTH > 0)
        cairo_rel_line_to(cv->cr, -(width - EXON_ARROW_WIDTH), 0);
      cairo_line_to(cv->cr, x, y + height / 2);
      cairo_line_to(cv->cr, x + MIN(width, EXON_ARROW_WIDTH), y + height);
      if (width - EXON_ARROW_WIDTH > 0)
        cairo_line_to(cv->cr, x + width, y + height);
      cairo_close_path(cv->cr);
      break;
    case STRAND_BOTH:
    case STRAND_UNKNOWN:
      cairo_rectangle(cv->cr, x, y, width, height);
   }

   cairo_fill_preserve(cv->cr);
   cairo_set_source_rgb(cv->cr, 0, 0, 0);
   cairo_stroke(cv->cr);
}

static void show_children(Genome_node *gn, void *data)
{
  Genome_feature *gf = (Genome_feature*) gn;
  Range feature_range = genome_node_get_range(gn);
  Show_children_info *info = (Show_children_info*) data,
                     local_info;
  char buf[BUFSIZ];
  double x_left, x_right, width, text_y;

  assert(info->cv->drawed_sequence_range_is_defined);
  assert(info->local_track_number < info->cv->number_of_tracks);

  text_y = info->local_track_number * TRACK_HEIGHT + TEXT_POSITION;
  x_left = info->cv->width * SPACE + info->cv->width * ROOM *
           ((double) (feature_range.start -
                      info->cv->drawed_sequence_range.start + 1) /
            (double) range_length(info->cv->drawed_sequence_range));
  x_right = info->cv->width * SPACE + info->cv->width * ROOM *
            ((double) (feature_range.end -
                       info->cv->drawed_sequence_range.start + 1) /
             (double) range_length(info->cv->drawed_sequence_range));
  width = info->cv->width * ROOM *
          ((double) range_length(feature_range) /
           (double) range_length(info->cv->drawed_sequence_range));

  cairo_set_source_rgb(info->cv->cr, 0, 0, 0);

  cairo_move_to(info->cv->cr, x_left, text_y);
  (void) snprintf(buf, BUFSIZ, "%lu", feature_range.start);
  cairo_show_text(info->cv->cr, buf);

  cairo_move_to(info->cv->cr, x_left + width / 2, text_y);
  cairo_show_text(info->cv->cr,
                  genome_feature_type_get_cstr(genome_feature_get_type(gf)));

  cairo_move_to(info->cv->cr, x_right, text_y);
  (void) snprintf(buf, BUFSIZ, "%lu", feature_range.end);
  cairo_show_text(info->cv->cr, buf);

  draw_exon_box(info->cv, gf, width, info->local_track_number);

  if (info->children_overlap)
    info->local_track_number++;

  if (genome_node_has_children(gn)) {
    local_info.cv = info->cv;
    local_info.children_overlap =
      !genome_node_direct_children_do_not_overlap(gn);
    local_info.local_track_number = info->local_track_number;
    if (!info->children_overlap)
      local_info.local_track_number++;

    genome_node_traverse_direct_children(gn, &local_info, show_children);

    if (!local_info.children_overlap)
      local_info.local_track_number++;
    assert(local_info.local_track_number >= info->local_track_number);

    info->local_track_number = local_info.local_track_number;

  if (!info->children_overlap)
    info->local_track_number--;
  }
}

static void cairo_visitor_genome_feature(Genome_visitor *gv,
                                         Genome_feature *gf,
                                         /*@unused@*/ Log *l)
{
  Cairo_visitor *cv = cairo_visitor_cast(gv);
  Show_children_info info;
  Range feature_range;
  char buf[BUFSIZ];
  double x_left, x_right, width, text_y;

  assert(cv->drawed_sequence_range_is_defined);

  feature_range = genome_node_get_range((Genome_node*) gf);

  /* reset track number if necessary */
  if (cv->last_range_is_defined &&
      !range_overlap(cv->last_range, feature_range)) {
    cv->global_track_number = 1;
    cv->last_range_is_defined = 0;
  }

  assert(cv->global_track_number < cv->number_of_tracks);

  /* calculate coordinates */
  x_left = cv->width * SPACE +
           cv->width * ROOM *
           ((double) (feature_range.start -
                      cv->drawed_sequence_range.start + 1)/
            (double) range_length(cv->drawed_sequence_range));
  x_right = cv->width * SPACE +
            cv->width * ROOM *
            ((double) (feature_range.end -
                       cv->drawed_sequence_range.start + 1)/
             (double) range_length(cv->drawed_sequence_range)),
  width = cv->width * ROOM *
          ((double) range_length(feature_range) /
           (double) range_length(cv->drawed_sequence_range));
  text_y = cv->global_track_number * TRACK_HEIGHT + TEXT_POSITION;

  cairo_set_source_rgb(cv->cr, 0, 0, 0);
  /* show <start --- feature_type --- end> */
  cairo_move_to(cv->cr, x_left, text_y);
  (void) snprintf(buf, BUFSIZ, "%lu", genome_node_get_start((Genome_node*) gf));
  cairo_show_text(cv->cr, buf);

  cairo_move_to(cv->cr, x_left + width / 2, text_y);
  cairo_show_text(cv->cr,
                  genome_feature_type_get_cstr(genome_feature_get_type(gf)));

  cairo_move_to(cv->cr, x_right, text_y);
  (void) snprintf(buf, BUFSIZ, "%lu", genome_node_get_end((Genome_node*) gf));
  cairo_show_text(cv->cr, buf);
  /* draw feature line */
  cairo_move_to(cv->cr, x_left,
                cv->global_track_number * TRACK_HEIGHT + FEATURE_POSITION);
  cairo_rel_line_to(cv->cr, width, 0);
  cairo_stroke(cv->cr);

  cv->global_track_number++;

  if (genome_node_has_children((Genome_node*) gf)) {
    info.cv = cv;
    info.children_overlap =
      !genome_node_direct_children_do_not_overlap((Genome_node*) gf);
    info.local_track_number = cv->global_track_number;
    genome_node_traverse_direct_children((Genome_node*) gf, &info,
                                         show_children);
    if (!info.children_overlap)
      info.local_track_number++;
    assert(info.local_track_number >= cv->global_track_number);
    cv->global_track_number = info.local_track_number;
  }

  /* store range for later use */
  if (cv->last_range_is_defined) {
    assert(feature_range.start >= cv->last_range.start);
    if (feature_range.end > cv->last_range.end)
      cv->last_range.end = feature_range.end;
  }
  else {
    cv->last_range = feature_range;
  }
  cv->last_range_is_defined = 1;
}

static void cairo_visitor_sequence_region(Genome_visitor *gv,
                                          Sequence_region *sr,
                                          /*@unused@*/ Log *l)
{
  Cairo_visitor *cv = cairo_visitor_cast(gv);
  Range sr_range;
  char buf[BUFSIZ];

  assert(!cv->drawed_sequence_range_is_defined);

  sr_range = genome_node_get_range((Genome_node*) sr);
  cv->drawed_sequence_range.start = MAX(sr_range.start, cv->from);
  cv->drawed_sequence_range.end = MIN(sr_range.end, cv->to);
  cv->drawed_sequence_range_is_defined = 1;

  assert(!cv->global_track_number);
  cairo_set_source_rgb(cv->cr, 0, 0, 0);
  /* show <start --- sequence id --- end> */
  cairo_move_to(cv->cr, cv->width * SPACE,
                cv->global_track_number * TRACK_HEIGHT + TEXT_POSITION);
  (void) snprintf(buf, BUFSIZ, "%lu", cv->drawed_sequence_range.start);
  cairo_show_text(cv->cr, buf);
  cairo_move_to(cv->cr, cv->width * (SPACE + ROOM / 2),
                cv->global_track_number * TRACK_HEIGHT + TEXT_POSITION);
  cairo_show_text(cv->cr, str_get(genome_node_get_seqid((Genome_node*) sr)));
  cairo_move_to(cv->cr, cv->width * (SPACE + ROOM),
                cv->global_track_number * TRACK_HEIGHT + TEXT_POSITION);
  (void) snprintf(buf, BUFSIZ, "%lu", cv->drawed_sequence_range.end);
  cairo_show_text(cv->cr, buf);
  /* draw sequence line */
  cairo_move_to(cv->cr, cv->width * SPACE,
                cv->global_track_number * TRACK_HEIGHT + FEATURE_POSITION);
  cairo_rel_line_to(cv->cr, cv->width * ROOM, 0);
  cairo_stroke(cv->cr);

  cv->global_track_number++;
}

const Genome_visitor_class* cairo_visitor_class()
{
  static Genome_visitor_class gvc = { sizeof (Cairo_visitor),
                                      cairo_visitor_free,
                                      NULL,
                                      cairo_visitor_genome_feature,
                                      cairo_visitor_sequence_region,
                                      NULL };
  return &gvc;
}

Genome_visitor* cairo_visitor_new(char *png_filename, int width,
                                  unsigned int number_of_tracks,
                                  unsigned long from, unsigned long to)
{
  Genome_visitor *gv = genome_visitor_create(cairo_visitor_class());
  Cairo_visitor *cv = cairo_visitor_cast(gv);

  assert(png_filename && width != UNDEFINT);

  cv->number_of_tracks = number_of_tracks;
  cv->from = from;
  cv->to = to;
  cv->width = width ;
  cv->height = cv->number_of_tracks * TRACK_HEIGHT;
  cv->global_track_number = 0;
  cv->surf = cairo_image_surface_create(CAIRO_FORMAT_RGB24, cv->width,
                                        cv->height);
  cv->cr = cairo_create(cv->surf);
  assert(cairo_status(cv->cr) == CAIRO_STATUS_SUCCESS);
  cv->png_filename = png_filename;
  cv->drawed_sequence_range_is_defined = 0;
  cv->last_range_is_defined = 0;

  cairo_set_source_rgb(cv->cr, 1, 1, 1);
  cairo_set_operator(cv->cr, CAIRO_OPERATOR_SOURCE);
  cairo_paint(cv->cr);

  return gv;
}
