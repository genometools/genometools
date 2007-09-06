/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef GENOME_NODE_INFO_H
#define GENOME_NODE_INFO_H

#include <stdbool.h>

/* some information about genome nodes stored efficiently in a single byte */
typedef char GenomeNodeInfo;

typedef enum {
  GENOME_NODE_TREE_STATUS_UNDETERMINED,
  GENOME_NODE_IS_TREE,
  GENOME_NODE_IS_NOT_A_TREE
} GenomeNodeTreeStatus;

/* add a parent */
void                 genome_node_info_add_parent(GenomeNodeInfo*);
/* return true if the genome node has multiple parents, false otherwise */
bool                 genome_node_info_multiple_parents(GenomeNodeInfo*);

/* return the tree status of the genome node */
GenomeNodeTreeStatus genome_node_info_get_tree_status(GenomeNodeInfo*);

void                 genome_node_info_set_tree_status(GenomeNodeInfo*,
                                                      GenomeNodeTreeStatus);

#endif
