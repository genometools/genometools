/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef MULTISET_MATCHING_H
#define MULTISET_MATCHING_H

/* Find all positions where a multiset defined by <multiset_string> occurs in
   text. For each such position the <ProcMatch> function is called.
   Positions start at 0. */
void multiset_matching(unsigned char *multiset_string,
                       unsigned long multiset_size, unsigned char *text,
                       unsigned long text_length, void *data,
                       void (*ProcMatch)(unsigned long pos, void *data));

#endif
