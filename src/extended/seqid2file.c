/*
  Copyright (c) 2007-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/seqid2file.h"

struct GtSeqid2FileInfo {
  GtStrArray *seqfiles;
  bool matchdesc,
       usedesc;
  GtStr *seqfile,
        *region_mapping;
};

GtSeqid2FileInfo* gt_seqid2file_info_new(void)
{
  GtSeqid2FileInfo *s2fi = gt_calloc(1, sizeof *s2fi);
  s2fi->seqfiles = gt_str_array_new();
  s2fi->seqfile = gt_str_new();
  s2fi->region_mapping = gt_str_new();
  return s2fi;
}

void gt_seqid2file_info_delete(GtSeqid2FileInfo *s2fi)
{
  if (!s2fi) return;
  gt_str_delete(s2fi->region_mapping);
  gt_str_delete(s2fi->seqfile);
  gt_str_array_delete(s2fi->seqfiles);
  gt_free(s2fi);
}

static int seqid2file_check(void *data, GT_UNUSED GtError *err)
{
  GtSeqid2FileInfo *info = (GtSeqid2FileInfo*) data;
  gt_error_check(err);
  gt_assert(info);
  if (gt_str_length(info->seqfile)) {
    gt_assert(!gt_str_array_size(info->seqfiles));
    gt_str_array_add(info->seqfiles, info->seqfile);
  }
  return 0;
}

void gt_seqid2file_register_options(GtOptionParser *op, GtSeqid2FileInfo *s2fi)
{
  GtOption *seqfile_option, *seqfiles_option, *matchdesc_option,
           *usedesc_option, *region_mapping_option;
  gt_assert(op && s2fi);

  /* -seqfile */
  seqfile_option = gt_option_new_filename("seqfile", "set the sequence file "
                                          "from which to extract the features",
                                          s2fi->seqfile);
  gt_option_parser_add_option(op, seqfile_option);

  /* -seqfiles */
  seqfiles_option = gt_option_new_filename_array("seqfiles", "set the sequence "
                                                 "files from which to extract "
                                                 "the features\nuse '--' to "
                                                 "terminate the list of "
                                                 "sequence files ",
                                                 s2fi->seqfiles);
  gt_option_parser_add_option(op, seqfiles_option);

  /* -matchdesc */
  matchdesc_option = gt_option_new_bool("matchdesc", "match the sequence "
                                        "descriptions from the input files for "
                                        "the desired sequence IDs (in GFF3)",
                                        &s2fi->matchdesc, false);
  gt_option_parser_add_option(op, matchdesc_option);

  /* -usedesc */
  usedesc_option = gt_option_new_bool("usedesc", "use sequence descriptions to "
                                      "map the sequence IDs (in GFF3) to "
                                      "actual sequence entries.\nIf a "
                                      "description contains a sequence range "
                                      "(e.g., III:1000001..2000000), the first "
                                      " part is used as sequence ID ('III') "
                                      "and the first range position as offset "
                                      "('1000001')", &s2fi->usedesc, false);
  gt_option_parser_add_option(op, usedesc_option);

  /* -regionmapping */
  region_mapping_option = gt_option_new_string("regionmapping", "set file "
                                               "containing sequence-region to "
                                               "sequence file mapping",
                                               s2fi->region_mapping, NULL);
  gt_option_parser_add_option(op, region_mapping_option);

  /* either option -seqfile, -seqfiles or -regionmapping is mandatory */
  gt_option_is_mandatory_either_3(seqfile_option, seqfiles_option,
                                  region_mapping_option);

  /* the options -seqfile and -regionmapping exclude each other */
  gt_option_exclude(seqfile_option, region_mapping_option);

  /* the options -seqfiles and -regionmapping exclude each other */
  gt_option_exclude(seqfiles_option, region_mapping_option);

  /* the options -seqfile and -seqfiles */
  gt_option_exclude(seqfile_option, seqfiles_option);

  /* the options -matchdesc and -usedesc exclude each other */
  gt_option_exclude(matchdesc_option, usedesc_option);

  /* option -matchdesc implies option -seqfile or -seqfiles*/
  gt_option_imply_either_2(matchdesc_option, seqfile_option, seqfiles_option);

  /* option -usedesc implies option -seqfile or -seqfiles */
  gt_option_imply_either_2(usedesc_option, seqfile_option, seqfiles_option);

  /* set hook function */
  gt_option_parser_register_hook(op, seqid2file_check, s2fi);
}

GtRegionMapping* gt_seqid2file_region_mapping_new(GtSeqid2FileInfo *s2fi,
                                                  GtError *err)
{
  gt_error_check(err);
  gt_assert(s2fi);
  /* create region mapping */
  if (gt_str_array_size(s2fi->seqfiles)) {
    return gt_region_mapping_new_seqfiles(s2fi->seqfiles, s2fi->matchdesc,
                                          s2fi->usedesc);
  }
  else
    return gt_region_mapping_new_mapping(s2fi->region_mapping, err);
}
