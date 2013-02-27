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

#ifndef GENOME_NODE_REP_H
#define GENOME_NODE_REP_H

#include <stdio.h>
#include "core/dlist.h"
#include "core/hashmap.h"
#include "core/thread_api.h"
#include "extended/genome_node.h"

typedef void    (*GtGenomeNodeFreeFunc)(GtGenomeNode*);
typedef GtStr*  (*GtGenomeNodeSetSeqidFunc)(GtGenomeNode*);
/* Used to sort nodes. */
typedef GtStr*  (*GtGenomeNodeGetIdstrFunc)(GtGenomeNode*);
typedef GtRange (*GtGenomeNodeGetRangeFunc)(GtGenomeNode*);
typedef void    (*GtGenomeNodeSetRangeFunc)(GtGenomeNode*, const GtRange*);
typedef void    (*GtGenomeNodeChangeSeqidFunc)(GtGenomeNode*, GtStr*);
typedef int     (*GtGenomeNodeAcceptFunc)(GtGenomeNode*, GtNodeVisitor*,
                                          GtError*);

/* the ``genome node'' interface */
struct GtGenomeNodeClass
{
  size_t size;
  GtGenomeNodeFreeFunc free;
  GtGenomeNodeSetSeqidFunc get_seqid;
  GtGenomeNodeGetIdstrFunc get_idstr;
  GtGenomeNodeGetRangeFunc get_range;
  GtGenomeNodeSetRangeFunc set_range;
  GtGenomeNodeChangeSeqidFunc change_seqid;
  GtGenomeNodeAcceptFunc accept;
};

struct GtGenomeNode
{
  const GtGenomeNodeClass *c_class;
  GtStr *filename;
  GtHashmap *userdata; /* created on demand */
  /* GtGenomeNodes are very space critical, therefore we can justify a bit
     ifdef-hell here... */
#ifdef GT_THREADS_ENABLED
  GtRWLock *lock;
#endif
  unsigned int line_number,
               reference_count,
               userdata_nof_items;
};

const GtGenomeNodeClass* gt_genome_node_class_new(size_t size,
                                       GtGenomeNodeFreeFunc free,
                                       GtGenomeNodeSetSeqidFunc get_seqid,
                                       GtGenomeNodeGetIdstrFunc get_idstr,
                                       GtGenomeNodeGetRangeFunc get_range,
                                       GtGenomeNodeSetRangeFunc set_range,
                                       GtGenomeNodeChangeSeqidFunc change_seqid,
                                       GtGenomeNodeAcceptFunc accept);
GtGenomeNode* gt_genome_node_create(const GtGenomeNodeClass*);

#endif
