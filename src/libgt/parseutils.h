/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef PARSEUTILS_H
#define PARSEUTILS_H

#include "error.h"
#include "range.h"
#include "phase.h"
#include "strand.h"

/* enforces that start <= end */
int    parse_range(Range*, const char *start, const char *end,
                   unsigned long line_number, const char *filename, Error*);

/* returns UNDEFDOUBLE if strcmp(score, ".") == 0 */
double parse_score(const char *score,
                   unsigned long line_number, const char *filename, Error*);

Strand parse_strand(const char *strand,
                    unsigned long line_number, const char *filename, Error*);

Phase  parse_phase(const char *phase,
                   unsigned long line_number, const char *filename, Error*);

int    parse_int(const char *integer,
                 unsigned long line_number, const char *filename, Error*);

#endif
