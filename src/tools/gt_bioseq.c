/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

typedef struct {
  bool recreate,
       showfasta,
       gc_content,
       stat;
  unsigned long showseqnum,
                width;
} Bioseq_arguments;

static OPrval parse_options(int *parsed_args, Bioseq_arguments *arguments,
                            int argc, char **argv, Env *env)
{
  Option *option, *option_showfasta, *option_showseqnum, *option_width,
         *option_stat;
  OptionParser *op;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("[option ...] sequence_file [...]",
                         "Construct the Biosequence files for the given "
                         "sequence_file(s) (if necessary).", env);

  /* -recreate */
  option = option_new_bool("recreate", "recreate Biosequence files, even if "
                           "they exist already", &arguments->recreate, false,
                           env);
  option_parser_add_option(op, option, env);

  /* -showfasta */
  option_showfasta = option_new_bool("showfasta", "show sequences on stdout "
                                     "(in fasta format)", &arguments->showfasta,
                                     false, env);
  option_parser_add_option(op, option_showfasta, env);

  /* -showseqnum */
  option_showseqnum = option_new_ulong_min("showseqnum", "show sequence with "
                                           "given number on stdout (in fasta "
                                           "format)", &arguments->showseqnum,
                                           UNDEFULONG, 1, env);
  option_parser_add_option(op, option_showseqnum, env);

  /* -gc-content */
  option = option_new_bool("gc-content", "show GC-content on stdout (for DNA "
                           "files)", &arguments->gc_content, false, env);
  option_parser_add_option(op, option, env);

  /* -stat */
  option_stat = option_new_bool("stat", "show sequence statistics",
                                &arguments->stat, false, env);
  option_parser_add_option(op, option_stat, env);

  /* -width */
  option_width = option_new_ulong("width", "set output width for showing of "
                                  "sequences (0 disables formatting)",
                                  &arguments->width, 0, env);
  option_parser_add_option(op, option_width, env);

  /* option implications */
  option_imply_either_2(option_width, option_showfasta, option_showseqnum, env);

  /* option exclusions */
  option_exclude(option_showfasta, option_stat, env);
  option_exclude(option_showfasta, option_showseqnum, env);
  option_exclude(option_showseqnum, option_stat, env);

  /* parse */
  oprval = option_parser_parse_min_args(op, parsed_args, argc, argv,
                                        versionfunc, 1, env);
  option_parser_delete(op, env);

  return oprval;
}

int gt_bioseq(int argc, char *argv[], Env *env)
{
  Bioseq_arguments arguments;
  Bioseq *bioseq;
  int parsed_args, has_err = 0;
  env_error_check(env);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  assert(parsed_args < argc);

  /* option -showseqnum makes only sense if we got a single sequence file */
  if (arguments.showseqnum != UNDEFULONG && parsed_args + 1 != argc) {
    env_error_set(env, "option '-showseqnum' makes only sense with a single "
                   "sequence_file");
    has_err = -1;
  }

  while (!has_err && parsed_args < argc) {
    /* bioseq construction */
    if (arguments.recreate)
      bioseq = bioseq_new_recreate(argv[parsed_args], env);
    else
      bioseq = bioseq_new(argv[parsed_args], env);
    if (!bioseq)
      has_err = -1;

    /* output */
    if (!has_err && arguments.showfasta)
      bioseq_show_as_fasta(bioseq, arguments.width);

    if (!has_err && arguments.showseqnum != UNDEFULONG) {
      if (arguments.showseqnum > bioseq_number_of_sequences(bioseq)) {
        env_error_set(env, "argument '%lu' to option '-showseqnum' is too "
                      "large. The Biosequence contains only '%lu' sequences.",
                      arguments.showseqnum, bioseq_number_of_sequences(bioseq));
        has_err = -1;
      }
      if (!has_err) {
        bioseq_show_sequence_as_fasta(bioseq, arguments.showseqnum - 1,
                                      arguments.width);
      }
    }

    if (!has_err && arguments.gc_content)
      bioseq_show_gc_content(bioseq, env);

    if (!has_err && arguments.stat)
      bioseq_show_stat(bioseq);

    /* free */
    bioseq_delete(bioseq, env);

    parsed_args++;
  }

  return has_err;
}
