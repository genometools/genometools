/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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
#include "libgtext/chseqids_stream.h"
#include "libgtext/genome_stream_rep.h"
#include "libgtext/mapping.h"

struct ChseqidsStream {
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  Mapping *chseqids_mapping;
  Array *genome_node_buffer;
  unsigned long buffer_index;
  bool sequence_regions_processed;
};

#define chseqids_stream_cast(GS)\
        genome_stream_cast(chseqids_stream_class(), GS)

static int change_sequence_id(GenomeNode *gn, void *data, Env *env)
{
  Str *changed_seqid = data;
  env_error_check(env);
  assert(changed_seqid);
  genome_node_set_seqid(gn, changed_seqid, env);
  return 0;
}

int chseqids_stream_next_tree(GenomeStream *gs, GenomeNode **gn, Env *env)
{
  ChseqidsStream *cs;
  GenomeNode *node, **gn_a, **gn_b;
  Str *changed_seqid;
  unsigned long i;
  int rval, had_err = 0;
  env_error_check(env);
  cs = chseqids_stream_cast(gs);

  if (!cs->sequence_regions_processed) {
    while (!had_err) {
      if (!(had_err = genome_stream_next_tree(cs->in_stream, &node, env))) {
        if (node)
          array_add(cs->genome_node_buffer, node);
        else
          break;
        if (!(genome_node_cast(sequence_region_class(), node)))
          break; /* no more sequence regions */
      }
    }
    /* now the buffer contains only sequence regions (except the last entry)
       -> change sequence ids */
    for (i = 0; !had_err && i < array_size(cs->genome_node_buffer); i++) {
      node = *(GenomeNode**) array_get(cs->genome_node_buffer, i);
      if (genome_node_get_seqid(node)) {
        if  ((changed_seqid = mapping_map_string(cs->chseqids_mapping,
                                           str_get(genome_node_get_seqid(node)),
                                                 env))) {
          rval = genome_node_traverse_children(node, changed_seqid,
                                               change_sequence_id, true, env);
          assert(!rval); /* change_sequence_id() is sane */
          str_delete(changed_seqid);
        }
        else
          had_err = -1;
       }
    }
    /* sort them */
    if (!had_err)
      genome_nodes_sort(cs->genome_node_buffer);
    /* consolidate them */
    for (i = 1; !had_err && i + 1 < array_size(cs->genome_node_buffer); i++) {
      gn_a = array_get(cs->genome_node_buffer, i-1);
      gn_b = array_get(cs->genome_node_buffer, i);
      if (genome_nodes_are_equal_sequence_regions(*gn_a, *gn_b)) {
        sequence_regions_consolidate(*gn_b, *gn_a);
        genome_node_rec_delete(*gn_a);
        *gn_a = NULL;
      }
    }
    cs->sequence_regions_processed = true;
  }

  /* return non-null nodes from buffer */
  while (!had_err && cs->buffer_index < array_size(cs->genome_node_buffer)) {
    node = *(GenomeNode**) array_get(cs->genome_node_buffer, cs->buffer_index);
    cs->buffer_index++;
    if (node) {
      *gn = node;
      return had_err;
    }
  }

  if (!had_err)
    had_err = genome_stream_next_tree(cs->in_stream, gn, env);
  if (!had_err && *gn) {
    if (genome_node_get_seqid(*gn)) {
      changed_seqid = mapping_map_string(cs->chseqids_mapping,
                                         str_get(genome_node_get_seqid(*gn)),
                                         env);
      assert(changed_seqid); /* is always defined, because an undefined mapping
                                would be catched earlier */
      rval = genome_node_traverse_children(*gn, changed_seqid,
                                           change_sequence_id, true, env);
      assert(!rval); /* change_sequence_id() is sane */
      str_delete(changed_seqid);
    }
  }

  return had_err;
}

static void chseqids_stream_free(GenomeStream *gs)
{
  ChseqidsStream *cs;
  unsigned long i;
  cs = chseqids_stream_cast(gs);
  mapping_delete(cs->chseqids_mapping);
  for (i = cs->buffer_index; i < array_size(cs->genome_node_buffer); i++) {
    genome_node_rec_delete(*(GenomeNode**)
                           array_get(cs->genome_node_buffer, i));
  }
  array_delete(cs->genome_node_buffer);
}

const GenomeStreamClass* chseqids_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (ChseqidsStream),
                                         chseqids_stream_next_tree,
                                         chseqids_stream_free };
  return &gsc;
}

GenomeStream* chseqids_stream_new(GenomeStream *in_stream, Str *chseqids_file,
                                  Env *env)
{
  GenomeStream *gs;
  ChseqidsStream *cs;
  env_error_check(env);
  assert(in_stream && chseqids_file);
  assert(genome_stream_is_sorted(in_stream));
  gs = genome_stream_create(chseqids_stream_class(), false, env);
  cs = chseqids_stream_cast(gs);
  cs->in_stream = in_stream;
  cs->chseqids_mapping = mapping_new(chseqids_file, "chseqids",
                                     MAPPINGTYPE_STRING, env_error(env));
  if (!cs->chseqids_mapping) {
    genome_stream_delete(gs);
    return NULL;
  }
  cs->genome_node_buffer = array_new(sizeof (GenomeNode*));
  return gs;
}
