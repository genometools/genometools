/*
  Copyright (c) 2007 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef CONFIG_H
#define CONFIG_H

#include <libgtcore/str.h>
#include <libgtext/color.h>

typedef struct Config Config;

Config*      config_new(Env*, bool*);
void         config_delete(Config*, Env*);
void         config_load_file(Config*, Str*, Env*);
void         config_reload(Config*, Env*);
Color        config_get_color(Config*, const char*, Env*);
void         config_set_color(Config*, const char*, Color, Env*);
void         config_set_cstr(Config *cfg,
                             const char* section,
														 const char *key,
														 const char* str,
														 Env* env);
const char*  config_get_cstr(Config *cfg,
                             const char* section,
														 const char *key,
														 Env* env);
int          config_unit_test(Env*);

#endif
