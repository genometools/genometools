/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef GENOME_NODE_REP_H
#define GENOME_NODE_REP_H

#include <stdio.h>
#include "libgtcore/dlist.h"
#include "libgtext/genome_node.h"

/* the ``genome node'' interface */
struct GenomeNodeClass
{
  size_t size;
  void  (*free)(GenomeNode*);
  Str*  (*get_seqid)(GenomeNode*);
  Str*  (*get_idstr)(GenomeNode*);
  Range (*get_range)(GenomeNode*);
  void  (*set_range)(GenomeNode*, Range);
  void  (*set_seqid)(GenomeNode*, Str*);
  int   (*accept)(GenomeNode*, GenomeVisitor*, Error*);
};

struct GenomeNode
{
  const GenomeNodeClass *c_class;
  Str *filename;
  Dlist *children;
  unsigned int line_number,
               reference_count,
               bit_field; /* uses the first 5 bits:
                             0:   mark
                             1-2: parent status
                             3-4: tree status
                          */
};

#define PARENT_STATUS_OFFSET  1
#define PARENT_STATUS_MASK    0x3
#define TREE_STATUS_OFFSET    3
#define TREE_STATUS_MASK      0x3

GenomeNode* genome_node_create(const GenomeNodeClass*);

#endif
