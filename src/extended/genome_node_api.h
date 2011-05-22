/*
  Copyright (c) 2006-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef GENOME_NODE_API_H
#define GENOME_NODE_API_H

#include "core/fptr_api.h"
#include "core/range_api.h"
#include "core/str_api.h"

typedef struct GtGenomeNodeClass GtGenomeNodeClass;

/* The <GtGenomeNode> interface. The different implementation of the
   <GtGenomeNode> interface represent different parts of genome annotations (as
   they are usually found in GFF3 files). */
typedef struct GtGenomeNode GtGenomeNode;

/* Increase the reference count for <genome_node> and return it.
   <genome_node> cannot be <NULL>.*/
GtGenomeNode* gt_genome_node_ref(GtGenomeNode *genome_node);

/* Decrease the reference count for <genome_node> or delete it, if this was the
   last reference. */
void          gt_genome_node_delete(GtGenomeNode *genome_node);

/* Return the sequence ID of <genome_node>.
   Corresponds to column 1 of GFF3 feature lines. */
GtStr*        gt_genome_node_get_seqid(GtGenomeNode *genome_node);

/* Return the genomic range of of <genome_node>.
   Corresponds to columns 4 and 5 of GFF3 feature lines. */
GtRange       gt_genome_node_get_range(GtGenomeNode *genome_node);

/* Return the start of <genome_node>.
   Corresponds to column 4 of GFF3 feature lines. */
unsigned long gt_genome_node_get_start(GtGenomeNode *genome_node);

/* Return the end of <genome_node>.
   Corresponds to column 5 of GFF3 feature lines. */
unsigned long gt_genome_node_get_end(GtGenomeNode *genome_node);

/* Return the length of <genome_node>.
   Computed from column 4 and 5 of GFF3 feature lines. */
unsigned long gt_genome_node_get_length(GtGenomeNode *genome_node);

/* Return the filename the <genome_node> was read from.
   If the node did not originate from a file, an appropriate string is
   returned. */
const char*   gt_genome_node_get_filename(const GtGenomeNode* genome_node);

/* Return the line of the source file the <genome_node> was encountered on
   (if the node was read from a file) */
unsigned int  gt_genome_node_get_line_number(const GtGenomeNode*);

/* Set the genomic range of <genome_node> to given <range>. */
void          gt_genome_node_set_range(GtGenomeNode *genome_node,
                                       const GtRange *range);

/* Attaches a pointer to <data> to the <node> using a given string as <key>. */
void          gt_genome_node_add_user_data(GtGenomeNode *node,
                                           const char *key,
                                           void *data,
                                           GtFree free_func);

/* Returns the pointer attached to the node for a given <key>. */
void*         gt_genome_node_get_user_data(const GtGenomeNode*,
                                           const char *key);

/* Calls the destructor function associated with the user data attached under
   the <key> on the attached data. */
void          gt_genome_node_release_user_data(GtGenomeNode*, const char *key);

#endif
