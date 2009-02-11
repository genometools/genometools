/*
  Copyright (c) 2006-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef NODE_STREAM_API_H
#define NODE_STREAM_API_H

#include <stdbool.h>
#include "core/error_api.h"
#include "extended/genome_node_api.h"

typedef struct GtNodeStreamClass GtNodeStreamClass;

/* The <GtNodeStream> interface. */
typedef struct GtNodeStream GtNodeStream;

/* Increase the reference count for <node_stream> and return it. */
GtNodeStream* gt_node_stream_ref(GtNodeStream *node_stream);
/* Try to get the the next <GtGenomeNode> from <node_stream> and store it in
   <genome_node> (transfers ownership to <genome_node>).
   If no error occurs, 0 is returned and <genome_node> contains either the next
   <GtGenomeNode> or <NULL>, if the <node_stream> is exhausted.
   If an error occurs, -1 is returned and <err> is set accordingly (the status
   of <genome_node> is undefined, but no ownership transfer occured). */
int           gt_node_stream_next(GtNodeStream *node_stream,
                                  GtGenomeNode **genome_node,
                                  GtError *err);
/* Calls <gt_node_stream_next()> on <node_stream> repeatedly until the
   <node_stream> is exhausted (0 is returned) or an error occurs (-1 is returned
   and <err> is set). All retrieved <GtGenomeNode>s are deleted automatically
   with calls to <gt_genome_node_delete()>.
   This method is basically a convenience method which simplifies calls to
   <gt_node_stream_next()> in a loop where retrieved <GtGenomeNode>s are not
   processed any further. */
int           gt_node_stream_pull(GtNodeStream *node_stream, GtError *err);
/* Return <true> if <node_stream> is a sorted stream, <false> otherwise. */
bool          gt_node_stream_is_sorted(GtNodeStream *node_stream);
/* Decrease the reference count for <node_stream> or delete it, if this was the
   last reference. */
void          gt_node_stream_delete(GtNodeStream *node_stream);

#endif
