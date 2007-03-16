/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef PARSEUTILS_H
#define PARSEUTILS_H

#include <libgtcore/env.h>
#include <libgtcore/range.h>
#include <libgtcore/phase.h>
#include <libgtcore/strand.h>

/* enforces that start <= end */
int parse_range(Range*, const char *start, const char *end,
                unsigned long line_number, const char *filename, Env*);

/* sets 'score_value' to UNDEFDOUBLE if strcmp(score, ".") == 0 */
int parse_score(double *score_value, const char *score,
                unsigned long line_number, const char *filename, Env*);

int parse_strand(Strand*, const char *strand,
                 unsigned long line_number, const char *filename, Env*);

int parse_phase(Phase*, const char *phase,
               unsigned long line_number, const char *filename, Env*);

int parse_int(int*, const char *integer,
              unsigned long line_number, const char *filename, Env*);

#endif
