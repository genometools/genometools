/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

typedef struct {
  bool recreate,
       showfasta,
       stat;
  unsigned long showseqnum,
                width;
} Bioseq_arguments;

static int parse_options(Bioseq_arguments *arguments, int argc, char **argv)
{
  Option *option, *option_showfasta, *option_showseqnum, *option_width,
         *option_stat;
  OptionParser *op;
  int parsed_args;
  op = option_parser_new("[option ...] sequence_file [...]",
                         "Construct the Biosequence files for the given "
                         "sequence_file(s) (if necessary).");

  /* -recreate */
  option = option_new_bool("recreate", "recreate Biosequence files, even if "
                           "they exist already", &arguments->recreate, false);
  option_parser_add_option(op, option);

  /* -showfasta */
  option_showfasta = option_new_bool("showfasta", "show sequences on stdout "
                                     "(in fasta format)", &arguments->showfasta,
                                     false);
  option_parser_add_option(op, option_showfasta);

  /* -showseqnum */
  option_showseqnum = option_new_ulong_min("showseqnum", "show sequence with "
                                           "given number on stdout (in fasta "
                                           "format)", &arguments->showseqnum,
                                           UNDEFULONG, 1);
  option_parser_add_option(op, option_showseqnum);

  /* -stat */
  option_stat = option_new_bool("stat", "show sequence statistics",
                                &arguments->stat, false);
  option_parser_add_option(op, option_stat);

  /* -width */
  option_width = option_new_ulong("width", "set output width for showing of "
                                  "sequences (0 disables formatting)",
                                  &arguments->width, 0);
  option_parser_add_option(op, option_width);

  /* option implications */
  option_imply_either_2(option_width, option_showfasta, option_showseqnum);

  /* option exclusions */
  option_exclude(option_showfasta, option_stat);
  option_exclude(option_showfasta, option_showseqnum);
  option_exclude(option_showseqnum, option_stat);

  /* parse */
  option_parser_parse_min_args(op, &parsed_args, argc, argv, versionfunc, 1);
  option_parser_free(op);

  return parsed_args;
}

int gt_bioseq(int argc, char *argv[])
{
  Bioseq_arguments arguments;
  Bioseq *bioseq;
  int parsed_args;

  /* option parsing */
  parsed_args = parse_options(&arguments, argc, argv);
  assert(parsed_args < argc);

  /* option -showseqnum makes only sense if we got a single sequence file */
  if (parsed_args + 1 != argc)
    error("option '-showseqnum' makes only sense with a single sequence_file");

  while (parsed_args < argc) {
    /* bioseq construction */
    bioseq = bioseq_new(argv[parsed_args]);
    bioseq_fill(bioseq, arguments.recreate);

    /* output */
    if (arguments.showfasta)
      bioseq_show_as_fasta(bioseq, arguments.width);

    if (arguments.showseqnum != UNDEFULONG) {
      if (arguments.showseqnum > bioseq_number_of_sequences(bioseq)) {
        error("argument '%lu' to option '-showseqnum' is too large. The "
              "Biosequence contains only '%lu' sequences.",
              arguments.showseqnum, bioseq_number_of_sequences(bioseq));
      }
      bioseq_show_sequence_as_fasta(bioseq, arguments.showseqnum - 1,
                                    arguments.width);
    }

    if (arguments.stat)
      bioseq_show_stat(bioseq);

    /* free */
    bioseq_free(bioseq);

    parsed_args++;
  }

  return EXIT_SUCCESS;
}
