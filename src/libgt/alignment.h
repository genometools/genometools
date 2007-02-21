/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <stdio.h>
#include "env.h"

/* the Alignment class (an Alignment object has to be contructed backwards) */
typedef struct Alignment Alignment;

Alignment*    alignment_new(void);
Alignment*    alignment_new_with_seqs(const char *u, unsigned long ulen,
                                      const char *v, unsigned long vlen);
void          alignment_set_seqs(Alignment*, const char *u, unsigned long ulen,
                                 const char *v, unsigned long vlen);
void          alignment_add_replacement(Alignment*);
void          alignment_add_deletion(Alignment*);
void          alignment_add_insertion(Alignment*);
void          alignment_remove_last(Alignment*); /* undo last add opteration */
unsigned long alignment_eval(const Alignment*); /* returns unit cost */
void          alignment_show(const Alignment*, FILE*);
void          alignment_show_multieop_list(const Alignment*, FILE*);
int           alignment_unit_test(Env*);
void          alignment_delete(Alignment*);

#endif
