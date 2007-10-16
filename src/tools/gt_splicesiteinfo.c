/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtcore/warning.h"
#include "libgtext/addintrons_stream.h"
#include "libgtext/gff3_in_stream.h"
#include "libgtext/splicesiteinfo_stream.h"
#include "libgtext/gtdatahelp.h"
#include "libgtext/seqid2file.h"

typedef struct {
  Str *seqfile,
      *regionmapping;
  bool addintrons;
} SpliceSiteInfoArguments;

static OPrval parse_options(int *parsed_args,
                            SpliceSiteInfoArguments *arguments, int argc,
                            const char **argv, Env *env)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("[option ...] [GFF3_file ...]", "Show information "
                         "about splice sites given in GFF3 files.", env);

  /* -seqfile and -regionmapping */
  seqid2file_options(op, arguments->seqfile, arguments->regionmapping, env);

  /* -addintrons */
  option = option_new_bool("addintrons", "add intron features between existing "
                           "exon features\n(before computing the information "
                           "to be shown)", &arguments->addintrons, false, env);
  option_parser_add_option(op, option, env);

  /* parse */
  option_parser_set_comment_func(op, gtdata_show_help, NULL);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, env);
  option_parser_delete(op, env);
  return oprval;
}

int gt_splicesiteinfo(int argc, const char **argv, Env *env)
{
  GenomeStream *gff3_in_stream,
               *addintrons_stream = NULL,
               *splicesiteinfo_stream = NULL;
  GenomeNode *gn;
  SpliceSiteInfoArguments arguments;
  RegionMapping *regionmapping;
  int parsed_args, had_err = 0;
  env_error_check(env);

  /* option parsing */
  arguments.seqfile = str_new(env);
  arguments.regionmapping = str_new(env);
  switch (parse_options(&parsed_args, &arguments, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
      str_delete(arguments.regionmapping, env);
      str_delete(arguments.seqfile, env);
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      str_delete(arguments.regionmapping, env);
      str_delete(arguments.seqfile, env);
      return 0;
  }

  if (!had_err) {
    /* create gff3 input stream */
    gff3_in_stream = gff3_in_stream_new_unsorted(argc - parsed_args,
                                                 argv + parsed_args, false,
                                                 false, env);

    /* create region mapping */
    regionmapping = seqid2file_regionmapping_new(arguments.seqfile,
                                                 arguments.regionmapping, env);
    if (!regionmapping)
      had_err = -1;
  }

  if (!had_err) {
    /* create addintrons stream (if necessary) */
    if (arguments.addintrons)
      addintrons_stream = addintrons_stream_new(gff3_in_stream, env);

    /* create extract feature stream */
    splicesiteinfo_stream = splicesiteinfo_stream_new(arguments.addintrons
                                                      ? addintrons_stream
                                                      : gff3_in_stream,
                                                      regionmapping, env);

    /* pull the features through the stream and free them afterwards */
    while (!(had_err = genome_stream_next_tree(splicesiteinfo_stream, &gn,
                                               env)) && gn) {
      genome_node_rec_delete(gn, env);
    }
  }

  if (!had_err) {
    if (!splicesiteinfo_stream_show(splicesiteinfo_stream, env)) {
      warning("input file(s) contained no intron, use option -addintrons to "
              "add introns automatically");
    }
  }

  /* free */
  genome_stream_delete(splicesiteinfo_stream, env);
  genome_stream_delete(addintrons_stream, env);
  genome_stream_delete(gff3_in_stream, env);
  str_delete(arguments.regionmapping, env);
  str_delete(arguments.seqfile, env);

  return had_err;
}
