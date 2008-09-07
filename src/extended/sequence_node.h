/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#ifndef SEQUENCE_NODE_H
#define SEQUENCE_NODE_H

/* implements the ``genome node'' interface */
typedef struct SequenceNode SequenceNode;

#include "core/str.h"
#include "extended/genome_node.h"

const GT_GenomeNodeClass* sequence_node_class(void);
GT_GenomeNode*            sequence_node_new(const char *description,
                                         Str *sequence); /* takes ownership */
const char*            sequence_node_get_description(const SequenceNode*);
const char*            sequence_node_get_sequence(const SequenceNode*);
unsigned long          sequence_node_get_sequence_length(const SequenceNode*);

#endif
