/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef OPTION_H
#define OPTION_H

#include <stdbool.h>
#include "libgtcore/range.h"
#include "libgtcore/str.h"
#include "libgtcore/strarray.h"

#define TERMINAL_WIDTH  80

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
typedef int  (*ShowCommentFunc)(const char *progname, void *data, Error*);
typedef int  (*OptionParserHookFunc)(void *data, Error*);

/* the option parser */
OptionParser* option_parser_new(const char *synopsis, const char *one_liner);
/* takes ownership */
void          option_parser_add_option(OptionParser*, Option*);
void          option_parser_set_comment_func(OptionParser*, ShowCommentFunc,
                                             void* data);
/* set the mailadress used in the final ``Report bugs to'' line of the -help
   output to <address>. It should be of the form "<bill@microsoft.com>" */
void          option_parser_set_mailaddress(OptionParser*, const char *address);
/* register a hook function. All registered hook functions are called at the end
   of option_parser_parse().
   This allows to have a module which registers a bunch of options in the option
   parser and automatically performs necessary postprocessing after the option
   parsing has been done via a hook function (see outputfile.[ch] for an
   example). */
void          option_parser_register_hook(OptionParser*, OptionParserHookFunc,
                                          void *data);
void          option_parser_set_min_args(OptionParser*, unsigned int);
void          option_parser_set_max_args(OptionParser*, unsigned int);
void          option_parser_set_min_max_args(OptionParser*, unsigned int,
                                                            unsigned int);
OPrval        option_parser_parse(OptionParser*, int *parsed_args, int argc,
                                  const char **argv, ShowVersionFunc, Error*);
void          option_parser_delete(OptionParser*);

/* the options

   option descriptions are automatically formatted to TERMINAL_WIDTH, but it is
   possible to embed newlines into the descriptions to manually affect the
   formating */
Option*        option_new_outputfile(FILE**);
Option*        option_new_verbose(bool *value);
Option*        option_new_debug(bool *value);
Option*        option_new_bool(const char *option_str, const char *description,
                               bool *value, bool default_value);
Option*        option_new_double(const char *option_str,
                                 const char *description, double *value,
                                 double default_value);
Option*        option_new_double_min(const char *option_str,
                                     const char *description, double *value,
                                     double default_value, double min_value);
Option*        option_new_double_min_max(const char *option_str,
                                         const char *description, double *value,
                                         double default_value, double min_value,
                                         double max_value);
/* enforces that the given argument is >= 0.0 and <= 1.0 */
Option*        option_new_probability(const char *option_str,
                                      const char *description , double *value,
                                      double default_value);
Option*        option_new_int(const char *option_str,
                              const char *description,
                              int *value, int default_value);
Option*        option_new_int_min(const char *option_str,
                                  const char *description, int *value,
                                  int default_value, int min_value);
Option*        option_new_int_max(const char *option_str,
                                  const char *description, int *value,
                                  int default_value, int max_value);
Option*        option_new_int_min_max(const char *option_str,
                                      const char *description,
                                      int *value, int default_value,
                                      int min_value, int max_value);
Option*        option_new_uint(const char *option_str, const char *description,
                               unsigned int *value, unsigned int default_value);
Option*        option_new_uint_min(const char *option_str,
                                   const char *description,
                                   unsigned int *value,
                                   unsigned int default_value,
                                   unsigned int min_value);
Option*        option_new_uint_max(const char *option_str,
                                   const char *description,
                                   unsigned int *value,
                                   unsigned int default_value,
                                   unsigned int max_value);
Option*        option_new_uint_min_max(const char *option_str,
                                       const char *description,
                                       unsigned int *value,
                                       unsigned int default_value,
                                       unsigned int min_value,
                                       unsigned int max_value);
Option*        option_new_long(const char *option_str, const char *description,
                               long *value, long default_value);
Option*        option_new_ulong(const char *option_str, const char *description,
                                unsigned long *value,
                                unsigned long default_value);
Option*        option_new_ulong_min(const char *option_str,
                                    const char *description,
                                    unsigned long *value,
                                    unsigned long default_value,
                                    unsigned long min_value);
Option*        option_new_ulong_min_max(const char *option_str,
                                        const char *description,
                                        unsigned long *value,
                                        unsigned long default_value,
                                        unsigned long min_value,
                                        unsigned long max_value);
/* if <default_value> equals NULL, UNDEF_LONG will be used as default for
   range->start and range->end */
Option*        option_new_range(const char *option_str, const char *description,
                                Range *value, Range *default_value);
Option*        option_new_range_min_max(const char *option_str,
                                        const char *description, Range *value,
                                        Range *default_value,
                                        unsigned long min_value,
                                        unsigned long max_value);
Option*        option_new_string(const char *option_str,
                                 const char *description,
                                 Str *value, const char *default_value);
Option*        option_new_stringarray(const char *option_str,
                                      const char *description, StrArray*);
/* add an option which allows only arguments given in the NULL terminated
   <domain> (<default_value> must be an entry of <domain> or NULL) */
Option*        option_new_choice(const char *option_str,
                                 const char *description, Str *value,
                                 const char *default_value,
                                 const char **domain);
Option*        option_new_filename(const char *option_str,
                                   const char *description, Str*);
Option*        option_new_filenamearray(const char *option_str,
                                        const char *description, StrArray*);
Option*        option_ref(Option*);
const char*    option_get_name(const Option *o);
void           option_is_mandatory(Option*);
void           option_is_mandatory_either(Option*, const Option*);
/* if this function was called, <o> is only shown in the output of -help+ */
void           option_is_extended_option(Option *o);
/* if this function was called, <o> is only shown in the output of -helpdev */
void           option_is_development_option(Option *o);
void           option_imply(Option*, const Option*);
void           option_imply_either_2(Option*, const Option*, const Option*);
void           option_exclude(Option*, Option*);
void           option_hide_default(Option*);
void           option_argument_is_optional(Option*);
bool           option_is_set(const Option*);
void           option_delete(Option*);

#endif
