/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <gtcore.h>
#include "cairo_stream.h"
#include "cairo_visitor.h"
#include "genome_stream_rep.h"

struct CairoStream {
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
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

static int determine_number_of_tracks(GenomeNode *gn, void *data, Env *env)
{
  unsigned int *number_of_tracks = (unsigned int*) data;

  env_error_check(env);

  if (genome_node_has_children(gn)) {
    if (genome_node_direct_children_do_not_overlap(gn, env)) {
      /* non-overlapping children need only one track */
      *number_of_tracks += 1;
    }
    else {
      /* overlapping children need a separate track for each one */
      *number_of_tracks += genome_node_number_of_children(gn);
    }

    /* recursion */
    genome_node_traverse_direct_children(gn, data, determine_number_of_tracks,
                                         env);
  }

  return 0;
}

static int cairo_stream_next_tree(GenomeStream *gs, GenomeNode **gn, Env *env)
{
  CairoStream *cairo_stream;
  GenomeNode *gn_ref;
  Range gn_range;
  int has_err;

  env_error_check(env);

  cairo_stream = cairo_stream_cast(gs);
  has_err = genome_stream_next_tree(cairo_stream->in_stream, gn, env);

  if (!has_err) {
    /* take the first one, if a seqid hasn't been defined already */
    if (!str_length(cairo_stream->seqid)) {
      str_delete(cairo_stream->seqid, env);
      cairo_stream->seqid = str_ref(genome_node_get_seqid(*gn));
    }

    assert(cairo_stream->seqid);
    /* XXX: we shown only features which lie completely in the given range */
    assert(cairo_stream->from);
    assert(cairo_stream->to);
    if ((str_cmp(cairo_stream->seqid, genome_node_get_seqid(*gn)) == 0) &&
        (!cairo_stream->sequence_region_added ||
         (genome_node_get_start(*gn) >= cairo_stream->from &&
          genome_node_get_end(*gn) <= cairo_stream->to))) {
      assert(cairo_stream->sequence_region_added ||
             genome_node_cast(sequence_region_class(), *gn));
      gn_ref = genome_node_rec_ref(*gn, env);
      array_add(cairo_stream->nodes_to_draw, gn_ref, env);

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
      has_err = determine_number_of_tracks(gn_ref,
                                           &cairo_stream->current_depth, env);
      if (!has_err) {
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
    }
  }

  assert(cairo_stream->number_of_tracks >= cairo_stream->current_depth);

  return has_err;
}

static void cairo_stream_free(GenomeStream *gs, Env *env)
{
  CairoStream *cairo_stream;
  unsigned long i;

  cairo_stream = cairo_stream_cast(gs);

  str_delete(cairo_stream->seqid, env);
  env_ma_free(cairo_stream->png_filename, env);
  for (i = 0; i < array_size(cairo_stream->nodes_to_draw); i++) {
   genome_node_rec_delete(*(GenomeNode**)
                          array_get(cairo_stream->nodes_to_draw, i), env);
  }
  array_delete(cairo_stream->nodes_to_draw, env);
}

const GenomeStreamClass* cairo_stream_class(void)
{
  static GenomeStreamClass gsc = { sizeof (CairoStream),
                                   cairo_stream_next_tree,
                                   cairo_stream_free };
  return &gsc;
}

GenomeStream* cairo_stream_new(GenomeStream *in_stream,
                                Str *seqid,
                                unsigned long from,
                                unsigned long to,
                                const char *png_filename,
                                int width, Env *env)
{
  CairoStream *cairo_stream;
  GenomeStream *gs;
  env_error_check(env);
  assert(seqid);
  assert(from <= to);
  gs = genome_stream_create(cairo_stream_class(),
                            genome_stream_is_sorted(in_stream), env);
  cairo_stream = cairo_stream_cast(gs);
  cairo_stream->in_stream = in_stream;
  cairo_stream->sequence_region_added = 0;
  cairo_stream->last_range_is_defined = 0;
  cairo_stream->seqid = str_ref(seqid);
  cairo_stream->from = from;
  cairo_stream->to = to;
  cairo_stream->current_depth = 0;
  cairo_stream->number_of_tracks = 0;
  cairo_stream->png_filename = cstr_dup(png_filename, env);
  cairo_stream->width = width;
  cairo_stream->nodes_to_draw = array_new(sizeof (GenomeNode*), env);
  return gs;
}

void cairo_stream_draw(CairoStream *cairo_stream, bool verbose, Env *env)
{
  GenomeVisitor *cairo_visitor;
  GenomeNode *gn;
  unsigned long i;

  env_error_check(env);
  assert(cairo_stream);

  cairo_visitor = cairo_visitor_new(cairo_stream->png_filename,
                                    cairo_stream->width,
                                    cairo_stream->number_of_tracks,
                                    cairo_stream->from,
                                    cairo_stream->to, env);

  /* init */
  i = 0;
  if (verbose) {
    printf("drawing file \"%s\"\n", cairo_stream->png_filename);
    progressbar_start(&i, array_size(cairo_stream->nodes_to_draw));
  }

  /* drawing */
  for (; i < array_size(cairo_stream->nodes_to_draw); i++) {
    gn = *(GenomeNode**) array_get(cairo_stream->nodes_to_draw, i);
    genome_node_accept(gn, cairo_visitor, env);
  }

  /* teardown */
  if (verbose)
    progressbar_stop();

  genome_visitor_delete(cairo_visitor, env);
}

