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
#include "cairo_stream.h"
#include "cairo_visitor.h"
#include "genome_stream_rep.h"
#include "progressbar.h"
#include "xansi.h"

struct _Cairo_stream {
  const Genome_stream parent_instance;
  Genome_stream *in_stream;
  unsigned int sequence_region_added : 1,
               last_range_is_defined : 1;
  Range last_range;
  Str *seqid; /* the id of the sequence region which is drawn */
  unsigned long from,
                to;
  char *png_filename;
  int width;
  unsigned int current_depth,
               number_of_tracks;
  Array *nodes_to_draw;
};

#define cairo_stream_cast(GS)\
        genome_stream_cast(cairo_stream_class(), GS);

static void determine_number_of_tracks(Genome_node *gn, void *data)
{
  unsigned int *number_of_tracks = (unsigned int*) data;

  if (genome_node_has_children(gn)) {
    if (genome_node_direct_children_do_not_overlap(gn)) {
      /* non-overlapping children need only one track */
      *number_of_tracks += 1;
    }
    else {
      /* overlapping children need a separate track for each one */
      *number_of_tracks += genome_node_number_of_children(gn);
    }

    /* recursion */
    genome_node_traverse_direct_children(gn, data, determine_number_of_tracks);
  }
}

static Genome_node* cairo_stream_next_tree(Genome_stream *gs, Log *l)
{
  Cairo_stream *cairo_stream;
  Genome_node *gn, *gn_ref;
  Range gn_range;

  cairo_stream = cairo_stream_cast(gs);
  gn = genome_stream_next_tree(cairo_stream->in_stream, l);

  if (!gn) return NULL;

  /* take the first one, if a seqid hasn't been defined already */
  if (!str_length(cairo_stream->seqid)) {
    str_free(cairo_stream->seqid);
    cairo_stream->seqid = str_ref(genome_node_get_seqid(gn));
  }

  assert(cairo_stream->seqid);
  /* XXX: we shown only features which lie completely in the given range */
  assert(cairo_stream->from);
  assert(cairo_stream->to);
  if ((str_cmp(cairo_stream->seqid, genome_node_get_seqid(gn)) == 0) &&
      (!cairo_stream->sequence_region_added ||
       (genome_node_get_start(gn) >= cairo_stream->from &&
        genome_node_get_end(gn) <= cairo_stream->to))) {
    assert(cairo_stream->sequence_region_added ||
           genome_node_cast(sequence_region_class(), gn));
    gn_ref = genome_node_rec_ref(gn);
    array_add(cairo_stream->nodes_to_draw, gn_ref);

    if (cairo_stream->last_range_is_defined &&
        range_overlap(cairo_stream->last_range,
                      genome_node_get_range(gn_ref))) {
      cairo_stream->current_depth++;
    }
    else {
      if (cairo_stream->sequence_region_added) {
        /* one for the sequence region and one for this feature */
        cairo_stream->current_depth = 2;
      }
      else {
        /* only this feature */
        cairo_stream->current_depth = 1;
      }
      cairo_stream->last_range_is_defined = 0;
    }
    determine_number_of_tracks(gn_ref, &cairo_stream->current_depth);
    if (cairo_stream->current_depth > cairo_stream->number_of_tracks)
      cairo_stream->number_of_tracks = cairo_stream->current_depth;

    gn_range = genome_node_get_range(gn_ref);
    if (cairo_stream->last_range_is_defined) {
      /* update range */
      assert(gn_range.start >= cairo_stream->last_range.start);
      if (gn_range.end > cairo_stream->last_range.end)
        cairo_stream->last_range.end = gn_range.end;
    }
    else if (cairo_stream->sequence_region_added) {
      /* save range */
      cairo_stream->last_range = gn_range;
      cairo_stream->last_range_is_defined = 1;
    }
    cairo_stream->sequence_region_added = 1;
  }

  assert(cairo_stream->number_of_tracks >= cairo_stream->current_depth);

  return gn;
}

static void cairo_stream_free(Genome_stream *gs)
{
  Cairo_stream *cairo_stream;
  unsigned long i;

  cairo_stream = cairo_stream_cast(gs);

  str_free(cairo_stream->seqid);
  free(cairo_stream->png_filename);
  for (i = 0; i < array_size(cairo_stream->nodes_to_draw); i++) {
   genome_node_rec_free(*(Genome_node**)
                        array_get(cairo_stream->nodes_to_draw, i));
  }
  array_free(cairo_stream->nodes_to_draw);
}

const Genome_stream_class* cairo_stream_class(void)
{
  static Genome_stream_class gsc = { sizeof (Cairo_stream),
                                     cairo_stream_next_tree,
                                     cairo_stream_free };
  return &gsc;
}

Genome_stream* cairo_stream_new(Genome_stream *in_stream,
                                Str *seqid,
                                unsigned long from,
                                unsigned long to,
                                const char *png_filename,
                                int width)
{
  Cairo_stream *cairo_stream;
  Genome_stream *gs = genome_stream_create(cairo_stream_class(),
                                           genome_stream_is_sorted(in_stream));
  assert(seqid);
  assert(from <= to);
  cairo_stream = cairo_stream_cast(gs);
  cairo_stream->in_stream = in_stream;
  cairo_stream->sequence_region_added = 0;
  cairo_stream->last_range_is_defined = 0;
  cairo_stream->seqid = str_ref(seqid);
  cairo_stream->from = from;
  cairo_stream->to = to;
  cairo_stream->current_depth = 0;
  cairo_stream->number_of_tracks = 0;
  cairo_stream->png_filename = xstrdup(png_filename);
  cairo_stream->width = width;
  cairo_stream->nodes_to_draw = array_new(sizeof (Genome_node*));
  return gs;
}

void cairo_stream_draw(Cairo_stream *cairo_stream, unsigned int verbose, Log *l)
{
  Genome_visitor *cairo_visitor;
  Genome_node *gn;
  unsigned long i;

  assert(cairo_stream);

  cairo_visitor = cairo_visitor_new(cairo_stream->png_filename,
                                    cairo_stream->width,
                                    cairo_stream->number_of_tracks,
                                    cairo_stream->from,
                                    cairo_stream->to);

  /* init */
  i = 0;
  if (verbose) {
    printf("drawing file \"%s\"\n", cairo_stream->png_filename);
    progressbar_start(&i, array_size(cairo_stream->nodes_to_draw));
  }

  /* drawing */
  for (; i < array_size(cairo_stream->nodes_to_draw); i++) {
    gn = *(Genome_node**) array_get(cairo_stream->nodes_to_draw, i);
    genome_node_accept(gn, cairo_visitor, l);
  }

  /* teardown */
  if (verbose)
    progressbar_stop();

  genome_visitor_free(cairo_visitor);
}

