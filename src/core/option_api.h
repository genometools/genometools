/*
  Copyright (c) 2006-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef OPTION_API_H
#define OPTION_API_H

#include <stdbool.h>
#include "core/range_api.h"
#include "core/str_api.h"
#include "core/str_array.h"

/* <GtOptionParser> objects can be used to parse command line options. */
typedef struct GtOptionParser GtOptionParser;
/* <GtOption> objects represent command line options (which are used in
   a <GtOptionParser>).
   Option descriptions are automatically formatted to
   <GT_OPTION_PARSER_TERMINAL_WIDTH>, but it is possible to embed newlines into
   the descriptions to manually affect the formating. */
typedef struct GtOption GtOption;

/* Possible option parser return values. <GT_OPTION_PARSER_OK> denotes that
   everything went fine, <GT_OPTION_PARSER_ERROR> that an error occured during
   option parsing, and <GT_OPTION_PARSER_REQUESTS_EXIT> that the option parser
   requests an exit, because option <-help>, <-help+>, <-helpdev> or <-version>
   was used. */
typedef enum GtOPrval GtOPrval;

enum GtOPrval{
  GT_OPTION_PARSER_OK,           /* Everything went fine. */
  GT_OPTION_PARSER_ERROR,        /* An error occured during option parsing. */
  GT_OPTION_PARSER_REQUESTS_EXIT /* The option parser requests an exit, because
                                    option -help, -help+, -helpdev, or -version
                                    was used. */
};

typedef void (*GtShowVersionFunc)(const char *progname);
typedef int  (*GtShowCommentFunc)(const char *progname, void *data, GtError*);
typedef int  (*GtOptionParserHookFunc)(void *data, GtError*);

/* The default terminal width used in the output of the <GtOptionParser>. */
#define GT_OPTION_PARSER_TERMINAL_WIDTH \
        80

/* Return a new <GtOptionParser> object. The <synopsis> should summarize the
   command line arguments and mandatory arguments in a single line.
   The <one_liner> should describe the program for which the <GtOptionParser> is
   used in a single line and must have an upper case letter at the start and a
   '.' at the end. */
GtOptionParser* gt_option_parser_new(const char *synopsis,
                                     const char *one_liner);
/* Add <option> to <option_parser>. Takes ownership of <option>. */
void            gt_option_parser_add_option(GtOptionParser *option_parser,
                                            GtOption *option);
/* Return the <GtOption> object if an option named <option_str> is present in
   <option_parser>, and <NULL> if no such option exists in <option_parser>. */
GtOption*       gt_option_parser_get_option(GtOptionParser *option_parser,
                                            const char *option_str);
/* Refer to manual at the end of <-help> output of <opion_parser>. */
void            gt_option_parser_refer_to_manual(GtOptionParser *option_parser);
/* Set <comment_func> in <option_parser> (<data> is passed along). */
void            gt_option_parser_set_comment_func(GtOptionParser *option_parser,
                                                  GtShowCommentFunc
                                                  comment_func,
                                                  void *data);
/* Set the version function used by <option_parser> to <version_func>.
   This version function takes precedence to the one supplied to
   <gt_option_parser_parse()>. */
void            gt_option_parser_set_version_func(GtOptionParser *option_parser,
                                                  GtShowVersionFunc
                                                  version_func);
/* Set the <mail_address> used in the final "Report bugs to" line of the <-help>
   output. It should be of the form <<bill@microsoft.com>> (email address
   enclosed in one pair of angle brackets). */
void            gt_option_parser_set_mail_address(GtOptionParser*,
                                                  const char *mail_address);
/* Register a <hook_function> with <option_parser>. All registered hook
   functions are called at the end of <gt_option_parser_parse(>).
   This allows to have a module which registers a bunch of options in the option
   parser and automatically performs necessary postprocessing after the option
   parsing has been done via a hook function. */
void            gt_option_parser_register_hook(GtOptionParser *option_parser,
                                               GtOptionParserHookFunc
                                               hook_function,
                                               void *data);
/* The the <minimum> number of additional command line arguments <option_parser>
   must parse in order to succeed. */
void            gt_option_parser_set_min_args(GtOptionParser *option_parser,
                                              unsigned int minimum);
/* The the <maximum> number of additional command line arguments <option_parser>
   must parse in order to succeed. */
void            gt_option_parser_set_max_args(GtOptionParser *option_parser,
                                              unsigned int maximum);
/* The the <minimum> and <maximum> number of additional command line arguments
   <option_parser> must parse in order to succeed. */
void            gt_option_parser_set_min_max_args(GtOptionParser *option_parser,
                                                  unsigned int minumum,
                                                  unsigned int maximum);
GtOPrval        gt_option_parser_parse(GtOptionParser *option_parser,
                                       int *parsed_args,
                                       int argc, const char **argv,
                                       GtShowVersionFunc version_func,
                                       GtError *err);
/* Delete <option_parser>. */
void            gt_option_parser_delete(GtOptionParser *option_parser);

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
/* Enforces that the given argument is >= 0.0 and <= 1.0. */
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
/* If <default_value> equals <NULL>, <GT_UNDEF_LONG> will be used as default for
   range->start and range->end. */
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
/* Add an option which allows only arguments given in the <NULL> terminated
   <domain> (<default_value> must be an entry of <domain> or <NULL>) */
GtOption*       gt_option_new_choice(const char *option_str,
                                     const char *description, GtStr *value,
                                     const char *default_value,
                                     const char **domain);
GtOption*       gt_option_new_filename(const char *option_str,
                                       const char *description, GtStr*);
GtOption*       gt_option_new_filenamearray(const char *option_str,
                                            const char *description,
                                            GtStrArray*);
/* Return a new debug <GtOption> object: <-debug>, "enable debugging output",
   default is <false>. The result of the option parsing is stored in <value> */
GtOption*       gt_option_new_debug(bool *value);
/* Return a new verbose <GtOption> object: <-v>, "be verbose",
   default is <false>. The result of the option parsing is stored in <value> */
GtOption*       gt_option_new_verbose(bool *value);
/* Return a new widht <GtOption> object: <-width>, "set output width for FASTA
   sequence printing (0 disables formatting)", default is 0.
   The result of the option parsing is stored in <value> */
GtOption*       gt_option_new_width(unsigned long *value);
/* Increate the reference count for <option> and return it. */
GtOption*       gt_option_ref(GtOption *option);
/* Return the name of <option> */
const char*     gt_option_get_name(const GtOption * option);
/* Make <option> mandatory. */
void            gt_option_is_mandatory(GtOption *option);
/* Make it mandatory, that either <option_a> or <option_b> is used. */
void            gt_option_is_mandatory_either(GtOption *option_a,
                                              const GtOption *option_b);
/* Make it mandatory, that one of the options <option_a>, <option_b>, or
   <option_c> is used. */
void            gt_option_is_mandatory_either_3(GtOption *option_a,
                                                const GtOption *option_b,
                                                const GtOption *option_c);
/* Set that <option> is only shown in the output of <-help+>. */
void            gt_option_is_extended_option(GtOption *option);
/* Set that <option> is only shown in the output of <-helpdev>. */
void            gt_option_is_development_option(GtOption *option);
/* Make <option_a> imply <option_b>. */
void            gt_option_imply(GtOption *option_a, const GtOption *option_b);
/* Make <option_b> imply either <option_b> or <option_c> */
void            gt_option_imply_either_2(GtOption *option_a,
                                         const GtOption *option_b,
                                         const GtOption *option_c);
/* Set that the options <option_a> and <option_b> exclude each other. */
void            gt_option_exclude(GtOption *option_a, GtOption *option_b);
/* Hide the default value of <option> in help output. */
void            gt_option_hide_default(GtOption *option);
/* Set that the argument to <option> is optional */
void            gt_option_argument_is_optional(GtOption *option);
/* Return <true> if <option> was set, <false> otherwise. */
bool            gt_option_is_set(const GtOption *option);
/* Delete <option>. */
void            gt_option_delete(GtOption*);

#endif
