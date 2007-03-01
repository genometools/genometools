/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "seqid2fileoptions.h"

void seqid2fileoptions(OptionParser *op, Str *seqfile, Str *regionmapping,
                       Env *env)
{
  Option *seqfile_option, *regionmapping_option;
  assert(op && seqfile && regionmapping);

  /* -seqfile */
  seqfile_option = option_new_string("seqfile", "set the sequence file from "
                                     "which to extract the features", seqfile,
                                     NULL, env);
  option_parser_add_option(op, seqfile_option, env);

  /* -regionmapping */
  regionmapping_option = option_new_string("regionmapping", "set file "
                                           "containing sequence-region to "
                                           "sequence file mapping",
                                           regionmapping, NULL, env);
  option_parser_add_option(op, regionmapping_option, env);

  /* either option -seqfile or -regionmapping is mandatory */
  option_is_mandatory_either(seqfile_option, regionmapping_option);

  /* the options -seqfile and -regionmapping exclude each other */
  option_exclude(seqfile_option, regionmapping_option, env);
}
