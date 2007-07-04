/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <gtcore.h>
#include <libgtview/png_stream.h>
#include <libgtview/png_visitor.h>
#include <libgtext/genome_stream_rep.h>

struct PNGStream {
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  bool sequence_region_added,
       last_range_is_defined;
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

#define png_stream_cast(GS)\
        genome_stream_cast(png_stream_class(), GS);

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

static int png_stream_next_tree(GenomeStream *gs, GenomeNode **gn, Env *env)
{
  PNGStream *png_stream;
  GenomeNode *gn_ref;
  Range gn_range;
  int had_err;

  env_error_check(env);

  png_stream = png_stream_cast(gs);
  had_err = genome_stream_next_tree(png_stream->in_stream, gn, env);

  if (!had_err && *gn) {
    /* take the first one, if a seqid hasn't been defined already */
    if (!str_length(png_stream->seqid)) {
      str_delete(png_stream->seqid, env);
      png_stream->seqid = str_ref(genome_node_get_seqid(*gn));
    }

    assert(png_stream->seqid);
    /* XXX: we shown only features which lie completely in the given range */
    assert(png_stream->from);
    assert(png_stream->to);
    if ((str_cmp(png_stream->seqid, genome_node_get_seqid(*gn)) == 0) &&
        (!png_stream->sequence_region_added ||
         (genome_node_get_start(*gn) >= png_stream->from &&
          genome_node_get_end(*gn) <= png_stream->to))) {
      assert(png_stream->sequence_region_added ||
             genome_node_cast(sequence_region_class(), *gn));
      gn_ref = genome_node_rec_ref(*gn, env);
      array_add(png_stream->nodes_to_draw, gn_ref, env);

      if (png_stream->last_range_is_defined &&
          range_overlap(png_stream->last_range,
                        genome_node_get_range(gn_ref))) {
        png_stream->current_depth++;
      }
      else {
        if (png_stream->sequence_region_added) {
          /* one for the sequence region and one for this feature */
          png_stream->current_depth = 2;
        }
        else {
          /* only this feature */
          png_stream->current_depth = 1;
        }
        png_stream->last_range_is_defined = false;
      }
      had_err = determine_number_of_tracks(gn_ref,
                                           &png_stream->current_depth, env);
      if (!had_err) {
        if (png_stream->current_depth > png_stream->number_of_tracks)
          png_stream->number_of_tracks = png_stream->current_depth;

        gn_range = genome_node_get_range(gn_ref);
        if (png_stream->last_range_is_defined) {
          /* update range */
          assert(gn_range.start >= png_stream->last_range.start);
          if (gn_range.end > png_stream->last_range.end)
            png_stream->last_range.end = gn_range.end;
        }
        else if (png_stream->sequence_region_added) {
          /* save range */
          png_stream->last_range = gn_range;
          png_stream->last_range_is_defined = true;
        }
        png_stream->sequence_region_added = true;
      }
    }
  }

  assert(png_stream->number_of_tracks >= png_stream->current_depth);

  return had_err;
}

static void png_stream_free(GenomeStream *gs, Env *env)
{
  PNGStream *png_stream;
  unsigned long i;

  png_stream = png_stream_cast(gs);

  str_delete(png_stream->seqid, env);
  env_ma_free(png_stream->png_filename, env);
  for (i = 0; i < array_size(png_stream->nodes_to_draw); i++) {
   genome_node_rec_delete(*(GenomeNode**)
                          array_get(png_stream->nodes_to_draw, i), env);
  }
  array_delete(png_stream->nodes_to_draw, env);
}

const GenomeStreamClass* png_stream_class(void)
{
  static GenomeStreamClass gsc = { sizeof (PNGStream),
                                   png_stream_next_tree,
                                   png_stream_free };
  return &gsc;
}

GenomeStream* png_stream_new(GenomeStream *in_stream, Str *seqid,
                             unsigned long from, unsigned long to,
                             const char *png_filename, int width, Env *env)
{
  PNGStream *png_stream;
  GenomeStream *gs;
  env_error_check(env);
  assert(seqid);
  assert(from <= to);
  gs = genome_stream_create(png_stream_class(),
                            genome_stream_is_sorted(in_stream), env);
  png_stream = png_stream_cast(gs);
  png_stream->in_stream = in_stream;
  png_stream->sequence_region_added = false;
  png_stream->last_range_is_defined = false;
  png_stream->seqid = str_ref(seqid);
  png_stream->from = from;
  png_stream->to = to;
  png_stream->current_depth = 0;
  png_stream->number_of_tracks = 0;
  png_stream->png_filename = cstr_dup(png_filename, env);
  png_stream->width = width;
  png_stream->nodes_to_draw = array_new(sizeof (GenomeNode*), env);
  return gs;
}

void png_stream_draw(PNGStream *png_stream, bool verbose, Env *env)
{
  GenomeVisitor *png_visitor;
  GenomeNode *gn;
  unsigned long i;

  env_error_check(env);
  assert(png_stream);

  png_visitor = png_visitor_new(png_stream->png_filename, png_stream->width,
                                png_stream->number_of_tracks, png_stream->from,
                                png_stream->to, env);

  /* init */
  i = 0;
  if (verbose) {
    printf("drawing file \"%s\"\n", png_stream->png_filename);
    progressbar_start(&i, array_size(png_stream->nodes_to_draw));
  }

  /* drawing */
  for (; i < array_size(png_stream->nodes_to_draw); i++) {
    gn = *(GenomeNode**) array_get(png_stream->nodes_to_draw, i);
    genome_node_accept(gn, png_visitor, env);
  }

  /* teardown */
  if (verbose)
    progressbar_stop();

  genome_visitor_delete(png_visitor, env);
}

