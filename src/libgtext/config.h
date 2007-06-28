/*
  Copyright (c) 2007 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef CONFIG_H
#define CONFIG_H

#include <libgtcore/str.h>
#include <libgtext/color.h>
#include <libgtext/genome_feature_type.h>

/* Represents domination status of an ordered pair.
   Used when two different types collapse into the same parent
   to determine splitting precedence. */
enum DominateStatus
{
  DOMINATES_FIRST,
  DOMINATES_SECOND,
  DOMINATES_EQUAL,
  DOMINATES_NOT_SPECIFIED,
  DOMINATES_UNKNOWN_TYPE
};

/* Holds configuration info for the gtview classes. */
typedef struct Config Config;

Config*      config_new(Env*, bool verbose);
void         config_load_file(Config*, Str*, Env*);
void         config_reload(Config*, Env*);
Color        config_get_color(Config*, const char*, Env*);
void         config_set_color(Config*, const char*, Color, Env*);
void         config_set_cstr(Config *cfg,
                             const char *section,
                             const char *key,
                             const char *str,
                             Env* env);
const char*  config_get_cstr(Config *cfg,
                             const char *section,
                             const char *key,
                             const char *deflt,
                             Env* env);
double       config_get_num(Config *cfg,
                            const char *section,
                            const char *key,
                            double deflt,
                            Env* env);
void         config_set_num(Config *cfg,
                            const char* section,
                            const char *key,
                            double number,
                            Env* env);
bool         config_cstr_in_list(Config *cfg,
                                 const char *section,
                                 const char *key,
                                 const char *checkstr,
                                 Env* env);
bool         config_get_verbose(Config *cfg);
int          config_dominates(Config *cfg,
                              GenomeFeatureType gft1,
                              GenomeFeatureType gft2,
                              Env* env);

int          config_unit_test(Env*);
void         config_delete(Config*, Env*);

#endif
