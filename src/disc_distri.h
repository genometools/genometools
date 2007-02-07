/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef DISC_DISTRI_H
#define DISC_DISTRI_H

/* A descrete distribution */
typedef struct Disc_distri Disc_distri;

Disc_distri* disc_distri_new(void);
void         disc_distri_add(Disc_distri*, unsigned long);
void         disc_distri_show(const Disc_distri*); /* on stdout */
void         disc_distri_free(Disc_distri*);

#endif
