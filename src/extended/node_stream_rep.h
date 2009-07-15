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

#ifndef NODE_STREAM_REP_H
#define NODE_STREAM_REP_H

#include <stdio.h>
#include "extended/node_stream_api.h"

typedef void (*GtNodeStreamFreeFunc)(GtNodeStream*);
typedef int  (*GtNodeStreamNextFunc)(GtNodeStream*, GtGenomeNode**, GtError*);

typedef struct GtNodeStreamMembers GtNodeStreamMembers;

struct GtNodeStream {
  const GtNodeStreamClass *c_class;
  GtNodeStreamMembers *members;
};

/* Create a new node stream class (that is, a class which implements the node
   stream interface). <size> denotes the size of objects of the new node stream
   class. The optional <free> method is called once, if an object of the new
   class is deleted. The mandatory <next> method has to implement the
   <gt_node_stream_next()> semantic for the new class. */
const
GtNodeStreamClass* gt_node_stream_class_new(size_t size,
                                            GtNodeStreamFreeFunc free,
                                            GtNodeStreamNextFunc next);

/* Create a new object of the given <node_stream_class>. If <ensure_sorting> is
   <true>, it is enforced that all genome node objects pulled from this class
   are sorted. That is, for consecutive nodes <a> and <b> obtained from the
   given <node_stream_class> the return code of <gt_genome_node_compare(a, b)>
   has to be smaller or equal than 0. If this condition is not met, an assertion
   fails. */
GtNodeStream*      gt_node_stream_create(const GtNodeStreamClass
                                         *node_stream_class,
                                         bool ensure_sorting);
/* Cast <node_stream> to the given <node_stream_class>.
   That is, if <node_stream> is not from the given <node_stream_class>, an
   assertion will fail. */
void*              gt_node_stream_cast(const GtNodeStreamClass
                                       *node_stream_class,
                                       GtNodeStream *node_stream);

#endif
