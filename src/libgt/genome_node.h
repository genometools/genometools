/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GENOME_NODE_H
#define GENOME_NODE_H

/* the ``genome node'' interface */
typedef struct Genome_node_class Genome_node_class;
typedef struct Genome_node Genome_node;

#include "bittab.h"
#include "genome_visitor.h"
#include "phase.h"
#include "range.h"
#include "str.h"

typedef void (*Genome_node_traverse_func)(Genome_node*, void*);

Genome_node*  genome_node_rec_ref(Genome_node*);
void*         genome_node_cast(const Genome_node_class*, Genome_node*);
void          genome_node_traverse_children(Genome_node*, void*,
                                            Genome_node_traverse_func,
                                            unsigned int traverse_only_once);
void          genome_node_traverse_direct_children(Genome_node*, void*,
                                                   Genome_node_traverse_func);
const char*   genome_node_get_filename(const Genome_node*);
unsigned long genome_node_get_line_number(const Genome_node*);
unsigned long genome_node_number_of_children(const Genome_node*);
Str*          genome_node_get_seqid(Genome_node*);
Str*          genome_node_get_idstr(Genome_node*); /* used to sort nodes */
unsigned long genome_node_get_start(Genome_node*);
unsigned long genome_node_get_end(Genome_node*);
Range         genome_node_get_range(Genome_node*);
void          genome_node_set_range(Genome_node*, Range);
void          genome_node_set_seqid(Genome_node*, Str*);
void          genome_node_set_source(Genome_node*, Str*);
void          genome_node_set_phase(Genome_node*, Phase);
void          genome_node_accept(Genome_node*, Genome_visitor*, Log *l);
void          genome_node_is_part_of_genome_node(Genome_node *parent,
                                                 Genome_node *child);
/* does not free the leaf */
void          genome_node_remove_leaf(Genome_node *tree, Genome_node *leafn);
unsigned int  genome_node_has_children(Genome_node*);
unsigned int  genome_node_direct_children_do_not_overlap(Genome_node*);
unsigned int  genome_node_tree_is_sorted(Genome_node **buffer,
                                         Genome_node *current_node);
/* returns 1 if the genome node overlaps at least one of the nodes given in the
   array. O(array_size) */
unsigned int  genome_node_overlaps_nodes(Genome_node*, Array*);
/* similar interface to genome_node_overlaps_nodes(). Aditionally, if a bittab
   is given (which must have the same size as the array), the bits corresponding
   to overlapped nodes are marked (i.e., set) */
unsigned int  genome_node_overlaps_nodes_mark(Genome_node*, Array*, Bittab*);
int           genome_node_compare(Genome_node**, Genome_node**);
void          genome_node_free(Genome_node*);
void          genome_node_rec_free(Genome_node*);

void          genome_nodes_sort(Array*);
void          genome_nodes_sort_stable(Array*);
unsigned int  genome_nodes_are_sorted(const Array*);

#endif
