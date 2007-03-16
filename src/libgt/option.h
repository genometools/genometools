/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef OPTION_H
#define OPTION_H

#include <stdbool.h>
#include <libgt/str.h>

/* the option parser class */
typedef struct OptionParser OptionParser;
/* the option class */
typedef struct Option Option;

/* possible option parser return values */
typedef enum {
  OPTIONPARSER_OK,           /* everything went fine */
  OPTIONPARSER_ERROR,        /* an error occured during option parsing */
  OPTIONPARSER_REQUESTS_EXIT /* the option parser requests an exit, because
                                option -help, -helpdev, or -version was used */
} OPrval;

typedef void (*ShowVersionFunc)(const char *progname);
typedef int  (*ShowCommentFunc)(const char *progname, void *data, Env*);

/* the option parser */
OptionParser* option_parser_new(const char *synopsis, const char *one_liner,
                                Env*);
/* takes ownership */
void          option_parser_add_option(OptionParser*, Option*, Env*);
void          option_parser_set_comment_func(OptionParser*, ShowCommentFunc,
                                             void* data);
OPrval        option_parser_parse(OptionParser*, int *parsed_args, int argc,
                                  const char **argv, ShowVersionFunc, Env*);
OPrval        option_parser_parse_min_args(OptionParser*, int *parsed_args,
                                           int argc, const char **argv,
                                           ShowVersionFunc,
                                           unsigned int
                                           min_additional_arguments, Env*);
OPrval        option_parser_parse_max_args(OptionParser*, int *parsed_args,
                                           int argc, const char **argv,
                                           ShowVersionFunc,
                                           unsigned int
                                           max_additional_arguments, Env*);
OPrval        option_parser_parse_min_max_args(OptionParser*, int *parsed_args,
                                               int argc, const char **argv,
                                               ShowVersionFunc,
                                               unsigned int
                                               min_additional_arguments,
                                               unsigned int
                                               max_additional_arguments,
                                               Env*);
void          option_parser_delete(OptionParser*, Env*);

/* the options */
Option*        option_new_outputfile(FILE**, Env*);
Option*        option_new_verbose(bool *value, Env*);
Option*        option_new_debug(bool *value, Env*);
Option*        option_new_bool(const char *option_str, const char *description,
                               bool *value, bool default_value, Env*);
Option*        option_new_double(const char *option_str,
                                 const char *description, double *value,
                                 double default_value, Env*);
Option*        option_new_int(const char *option_str,
                              const char *description,
                              int *value, int default_value, Env*);
Option*        option_new_int_min(const char *option_str,
                                  const char *description, int *value,
                                  int default_value, int min_value, Env*);
Option*        option_new_uint(const char *option_str, const char *description,
                               unsigned int *value, unsigned int default_value,
                               Env*);
Option*        option_new_uint_min(const char *option_str,
                                   const char *description,
                                   unsigned int *value,
                                   unsigned int default_value,
                                   unsigned int min_value, Env*);
Option*        option_new_uint_max(const char *option_str,
                                   const char *description,
                                   unsigned int *value,
                                   unsigned int default_value,
                                   unsigned int max_value, Env*);
Option*        option_new_long(const char *option_str, const char *description,
                               long *value, long default_value, Env*);
Option*        option_new_ulong(const char *option_str, const char *description,
                                unsigned long *value,
                                unsigned long default_value, Env*);
Option*        option_new_ulong_min(const char *option_str,
                                    const char *description,
                                    unsigned long *value,
                                    unsigned long default_value,
                                    unsigned long min_value, Env*);
Option*        option_new_string(const char *option_str,
                                 const char *description,
                                 Str *value, const char *default_value, Env*);
void           option_is_mandatory(Option*);
void           option_is_mandatory_either(Option*, const Option*);
void           option_is_development_option(Option*);
void           option_imply(Option*, const Option*, Env*);
void           option_imply_either_2(Option*, const Option*, const Option*,
                                     Env*);
void           option_exclude(Option*, Option*, Env*);
void           option_hide_default(Option*);
void           option_delete(Option*, Env*);

#endif
