/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef CSTR_H
#define CSTR_H

#include <stdio.h>
#include <libgtcore/env.h>
#include <libgtcore/genfile.h>

char*         cstr_dup(const char*, Env*);
void          cstr_show(const char*, unsigned long length, FILE*);

/* use p and a blank as prefix for cstr_array[0] and return the result */
char**        cstr_array_prefix_first(const char **cstr_array, const char *p,
                                      Env*);
char**        cstr_array_preprend(const char **cstr_array, const char *p, Env*);
void          cstr_array_show(char **cstr_array, FILE*);
void          cstr_array_show_genfile(const char **cstr_array, GenFile*);
unsigned long cstr_array_size(const char **cstr_array); /* O(n) */
void          cstr_array_delete(char **cstr_array, Env*);

#endif
