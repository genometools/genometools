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
#include "core/range.h"
#include "core/str.h"
#include "core/str_array.h"

#define TERMINAL_WIDTH  80

/* the option parser class */
typedef struct GtOptionParser GtOptionParser;
/* the option class */
typedef struct GtOption GtOption;

/* possible option parser return values */
typedef enum {
  OPTIONPARSER_OK,           /* everything went fine */
  OPTIONPARSER_ERROR,        /* an error occured during option parsing */
  OPTIONPARSER_REQUESTS_EXIT /* the option parser requests an exit, because
                                option -help, -helpdev, or -version was used */
} GtOPrval;

typedef void (*GtShowVersionFunc)(const char *progname);
typedef int  (*GtShowCommentFunc)(const char *progname, void *data, GtError*);
typedef int  (*GtOptionParserHookFunc)(void *data, GtError*);

/* the option parser */
GtOptionParser* gt_option_parser_new(const char *synopsis,
                                     const char *one_liner);
/* takes ownership */
void            gt_option_parser_add_option(GtOptionParser*, GtOption*);
/* refer to manual at the end of help output */
void            gt_option_parser_refer_to_manual(GtOptionParser*);
void            gt_option_parser_set_comment_func(GtOptionParser*,
                                                  GtShowCommentFunc,
                                                  void* data);
/* set the mailadress used in the final ``Report bugs to'' line of the -help
   output to <address>. It should be of the form "<bill@microsoft.com>" */
void            gt_option_parser_set_mailaddress(GtOptionParser*,
                                                 const char *address);
/* register a hook function. All registered hook functions are called at the end
   of gt_option_parser_parse().
   This allows to have a module which registers a bunch of options in the option
   parser and automatically performs necessary postprocessing after the option
   parsing has been done via a hook function (see outputfile.[ch] for an
   example). */
void            gt_option_parser_register_hook(GtOptionParser*,
                                               GtOptionParserHookFunc,
                                               void *data);
void            gt_option_parser_set_min_args(GtOptionParser*, unsigned int);
void            gt_option_parser_set_max_args(GtOptionParser*, unsigned int);
void            gt_option_parser_set_min_max_args(GtOptionParser*, unsigned int,
                                                  unsigned int);
GtOPrval        gt_option_parser_parse(GtOptionParser*, int *parsed_args,
                                       int argc, const char **argv,
                                       GtShowVersionFunc, GtError*);
void            gt_option_parser_delete(GtOptionParser*);

/* the options

   option descriptions are automatically formatted to TERMINAL_WIDTH, but it is
   possible to embed newlines into the descriptions to manually affect the
   formating */
GtOption*       gt_option_new_outputfile(FILE**);
GtOption*       gt_option_new_verbose(bool *value);
GtOption*       gt_option_new_debug(bool *value);
GtOption*       gt_option_new_bool(const char *option_str,
                                   const char *description,
                                   bool *value, bool default_value);
GtOption*       gt_option_new_double(const char *option_str,
                                     const char *description, double *value,
                                     double default_value);
GtOption*       gt_option_new_double_min(const char *option_str,
                                         const char *description, double *value,
                                         double default_value,
                                         double min_value);
GtOption*       gt_option_new_double_min_max(const char *option_str,
                                             const char *description,
                                             double *value,
                                             double default_value,
                                             double min_value,
                                             double max_value);
/* enforces that the given argument is >= 0.0 and <= 1.0 */
GtOption*       gt_option_new_probability(const char *option_str,
                                          const char *description,
                                          double *value,
                                          double default_value);
GtOption*       gt_option_new_int(const char *option_str,
                                  const char *description,
                                  int *value, int default_value);
GtOption*       gt_option_new_int_min(const char *option_str,
                                      const char *description, int *value,
                                      int default_value, int min_value);
GtOption*       gt_option_new_int_max(const char *option_str,
                                      const char *description, int *value,
                                      int default_value, int max_value);
GtOption*       gt_option_new_int_min_max(const char *option_str,
                                          const char *description,
                                          int *value, int default_value,
                                          int min_value, int max_value);
GtOption*       gt_option_new_uint(const char *option_str,
                                   const char *description,
                                   unsigned int *value,
                                   unsigned int default_value);
GtOption*       gt_option_new_uint_min(const char *option_str,
                                       const char *description,
                                       unsigned int *value,
                                       unsigned int default_value,
                                       unsigned int min_value);
GtOption*       gt_option_new_uint_max(const char *option_str,
                                       const char *description,
                                       unsigned int *value,
                                       unsigned int default_value,
                                       unsigned int max_value);
GtOption*       gt_option_new_uint_min_max(const char *option_str,
                                           const char *description,
                                           unsigned int *value,
                                           unsigned int default_value,
                                           unsigned int min_value,
                                           unsigned int max_value);
GtOption*       gt_option_new_long(const char *option_str,
                                   const char *description,
                                   long *value, long default_value);
GtOption*       gt_option_new_ulong(const char *option_str,
                                    const char *description,
                                    unsigned long *value,
                                    unsigned long default_value);
GtOption*       gt_option_new_ulong_min(const char *option_str,
                                        const char *description,
                                        unsigned long *value,
                                        unsigned long default_value,
                                        unsigned long min_value);
GtOption*       gt_option_new_ulong_min_max(const char *option_str,
                                            const char *description,
                                            unsigned long *value,
                                            unsigned long default_value,
                                            unsigned long min_value,
                                            unsigned long max_value);
/* if <default_value> equals NULL, GT_UNDEF_LONG will be used as default for
   range->start and range->end */
GtOption*       gt_option_new_range(const char *option_str,
                                    const char *description,
                                    GtRange *value, GtRange *default_value);
GtOption*       gt_option_new_range_min_max(const char *option_str,
                                            const char *description,
                                            GtRange *value,
                                            GtRange *default_value,
                                            unsigned long min_value,
                                            unsigned long max_value);
GtOption*       gt_option_new_string(const char *option_str,
                                     const char *description,
                                     GtStr *value, const char *default_value);
GtOption*       gt_option_new_stringarray(const char *option_str,
                                          const char *description, GtStrArray*);
/* add an option which allows only arguments given in the NULL terminated
   <domain> (<default_value> must be an entry of <domain> or NULL) */
GtOption*       gt_option_new_choice(const char *option_str,
                                     const char *description, GtStr *value,
                                     const char *default_value,
                                     const char **domain);
GtOption*       gt_option_new_filename(const char *option_str,
                                       const char *description, GtStr*);
GtOption*       gt_option_new_filenamearray(const char *option_str,
                                            const char *description,
                                            GtStrArray*);
GtOption*       gt_option_ref(GtOption*);
const char*     gt_option_get_name(const GtOption *o);
void            gt_option_is_mandatory(GtOption*);
void            gt_option_is_mandatory_either(GtOption*, const GtOption*);
/* if this function was called, <o> is only shown in the output of -help+ */
void            gt_option_is_extended_option(GtOption *o);
/* if this function was called, <o> is only shown in the output of -helpdev */
void            gt_option_is_development_option(GtOption *o);
void            gt_option_imply(GtOption*, const GtOption*);
void            gt_option_imply_either_2(GtOption*, const GtOption*,
                                         const GtOption*);
void            gt_option_exclude(GtOption*, GtOption*);
void            gt_option_hide_default(GtOption*);
void            gt_option_argument_is_optional(GtOption*);
bool            gt_option_is_set(const GtOption*);
void            gt_option_delete(GtOption*);

#endif
