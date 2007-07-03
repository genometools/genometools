/*
  Copyright (c) 2007 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef CONFIG_H
#define CONFIG_H

#include <libgtcore/str.h>
#include <libgtview/color.h>
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

/*!
Creates a Config object.
\param env Pointer to Environment object.
\param verbose Verbosity flag. If set, warnings will be given.
\return Pointer to the new object.
*/
Config*      config_new(Env *env, bool verbose);

/*!
Loads and executes a Lua configuration file.
This file must contain a global table called 'config'.
\param cfg Config object to load into.
\param fn Filename of the script to execute.
\param env Pointer to Environment object.
*/
void         config_load_file(Config *cfg, Str *fn, Env *env);

/*!
Reloads the Lua configuration file.
\param cfg Config object to search in.
\param env Pointer to Environment object.
*/
void         config_reload(Config *cfg, Env *env);

/*!
Retrieves a color value from the configuration.
\param cfg Config object to search in.
\param key Key (e.g. feature) to get a color for.
\param env Pointer to Environment object.
\return color Color associated with key.
*/
Color        config_get_color(Config *cfg,
                              const char *key,
                              Env *env);

/*!
Sets a color value in the configuration to a certain value.
\param cfg Config object to search in.
\param key Key (e.g. feature) to set a color for.
\param color Color to associate with key.
\param env Pointer to Environment object.
*/
void         config_set_color(Config *cfg,
                              const char *key,
                              Color color,
                              Env* env);

/*!
Sets a string value in the configuration to a certain value.
\param cfg Config object to search in.
\param section Section to set a key in.
\param key Key to set a value for.
\param str String that is to be set.
\param env Pointer to Environment object.
*/
void         config_set_cstr(Config *cfg,
                             const char *section,
                             const char *key,
                             const char *str,
                             Env* env);

/*!
Retrieves a string value from the configuration.
\param cfg Config object to search in.
\param section Section to get a key from.
\param key Key to get a value from.
\param deflt Default value to return if key not found.
\param env Pointer to Environment object.
\return string pointer to result, defaults to argument
*/
const char*  config_get_cstr(Config *cfg,
                             const char *section,
                             const char *key,
                             const char *deflt,
                             Env* env);

/*!
Retrieves a numeric value from the configuration.
\param cfg Config object to search in.
\param section Section to get a key from.
\param key Key to get a value from.
\param deflt Default value to return if key not found.
\param env Pointer to Environment object.
\return double result, defaults to argument
*/
double       config_get_num(Config *cfg,
                            const char *section,
                            const char *key,
                            double deflt,
                            Env* env);

/*!
Sets a numeric value in the configuration to a certain value.
\param cfg Config object to search in.
\param section Section to set a key in.
\param key Key to set a value for.
\param number Value that is to be set.
\param env Pointer to Environment object.
*/
void         config_set_num(Config *cfg,
                            const char* section,
                            const char *key,
                            double number,
                            Env* env);

/*!
Checks whether a given string appears in a list (table) of strings
in the configuration settings.
\param cfg Config object to search in.
\param cfg Section object to search in.
\param key Key (e.g. feature) to te list to be checked.
\param checkstr String whose membership is to be determined.
\param env Pointer to Environment object.
\return TRUE if checkstr is in list, FALSE otherwise
*/
bool         config_cstr_in_list(Config *cfg,
                                 const char *section,
                                 const char *key,
                                 const char *checkstr,
                                 Env* env);

/*!
Returns verbosity flag.
\param cfg Pointer to Config object.
\return Verbosity status as bool
*/bool         config_get_verbose(Config *cfg);

/*!
Compares two GenomeFeatureTypes w.r.t. their splitting
precendence as defined in the config object.
If a type dominates, it will be drawn on top of the other in the image.
\param cfg Pointer to Config object.
\param ft1 GFT to check
\param ft2 GFT to check
\param cfg Pointer to Environment object.
\return DominateStatus enum value
*/
int          config_dominates(Config *cfg,
                              GenomeFeatureType gft1,
                              GenomeFeatureType gft2,
                              Env* env);
/* Unit test */
int          config_unit_test(Env*);

/*!
Deletes a Config object.
\param cfg Pointer to Config object to delete.
\param env Pointer to Environment object.
*/
void         config_delete(Config*, Env*);

#endif
