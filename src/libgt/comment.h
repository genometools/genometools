/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef COMMENT_H
#define COMMENT_H

#include "str.h"

/* implements the ``genome node'' interface */
typedef struct Comment Comment;

#include "genome_node.h"

const GenomeNodeClass* comment_class(void);
GenomeNode*            comment_new(const char *comment, const char *filename,
                                   unsigned long line_number);
const char*            comment_get_comment(Comment *c);

#endif
