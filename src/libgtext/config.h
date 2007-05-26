/*
  Copyright (c) 2007 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef CONFIG_H
#define CONFIG_H

typedef struct Config Config;

typedef struct {
  double red, green, blue;
} Color;

Config*    config_new(Env*);
void       config_delete(Config*, Env*);
void       config_load_file(Config*, Str*, Env*);
void       config_reload(Config*, Env*);
Color      config_get_color(Config*, const char*, Env*);
void       config_set_color(Config*, const char*, Color, Env*);
int        config_unit_test(Env*);
#endif
