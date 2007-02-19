/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
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
