/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef SPLICEDSEQ_H
#define SPLICEDSEQ_H

typedef struct Splicedseq Splicedseq;

Splicedseq*   splicedseq_new(void);
/* adds an ``exon'' to the spliced sequence */
void          splicedseq_add(Splicedseq*, unsigned long start,
                             unsigned long end, const char *original_sequence);
char*         splicedseq_get(const Splicedseq*);
unsigned int  splicedseq_pos_is_border(const Splicedseq*, unsigned long);
/* maps the given position back to the original coordinate system */
unsigned long splicedseq_map(const Splicedseq*, unsigned long);
unsigned long splicedseq_length(const Splicedseq*);
void          splicedseq_reverse(Splicedseq*);
void          splicedseq_reset(Splicedseq*);
int           splicedseq_unit_test(void);
void          splicedseq_free(Splicedseq*);

#endif
