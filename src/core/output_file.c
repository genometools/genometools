/*
  Copyright (c) 2007-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <string.h>
#include "core/file.h"
#include "core/fileutils_api.h"
#include "core/ma.h"
#include "core/output_file_api.h"
#include "core/warning_api.h"

struct GtOutputFileInfo {
  GtStr *output_filename;
  bool gzip,
       bzip2,
       force;
  GtFile **outfp;
};

GtOutputFileInfo* gt_output_file_info_new(void)
{
  GtOutputFileInfo *ofi = gt_malloc(sizeof *ofi);
  ofi->output_filename = gt_str_new();
  return ofi;
}

static int determine_outfp(void *data, GtError *err)
{
  GtOutputFileInfo *ofi = (GtOutputFileInfo*) data;
  GtFileMode file_mode;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(ofi);
  if (!gt_str_length(ofi->output_filename))
    *ofi->outfp = NULL; /* no output file given -> use stdin */
  else { /* outputfile given -> create generic file pointer */
    gt_assert(!(ofi->gzip && ofi->bzip2));
    if (ofi->gzip)
      file_mode = GT_FILE_MODE_GZIP;
    else if (ofi->bzip2)
      file_mode = GT_FILE_MODE_BZIP2;
    else
      file_mode = GT_FILE_MODE_UNCOMPRESSED;
    if (file_mode != GT_FILE_MODE_UNCOMPRESSED &&
        strcmp(gt_str_get(ofi->output_filename) +
               gt_str_length(ofi->output_filename) -
               strlen(gt_file_mode_suffix(file_mode)),
               gt_file_mode_suffix(file_mode))) {
      gt_warning("output file '%s' doesn't have correct suffix '%s', appending "
                 "it", gt_str_get(ofi->output_filename),
                 gt_file_mode_suffix(file_mode));
      gt_str_append_cstr(ofi->output_filename, gt_file_mode_suffix(file_mode));
    }
    if (!ofi->force && gt_file_exists(gt_str_get(ofi->output_filename))) {
        gt_error_set(err, "file \"%s\" exists already, use option -%s to "
                     "overwrite", gt_str_get(ofi->output_filename),
                     GT_FORCE_OPT_CSTR);
        had_err = -1;
    }
    if (!had_err) {
      *ofi->outfp = gt_file_xopen_file_mode(file_mode,
                                            gt_str_get(ofi->output_filename),
                                            "w");
      gt_assert(*ofi->outfp);
    }
  }
  return had_err;
}

void gt_output_file_register_options(GtOptionParser *op, GtFile **outfp,
                                    GtOutputFileInfo *ofi)
{
  GtOption *opto, *optgzip, *optbzip2, *optforce;
  gt_assert(outfp && ofi);
  ofi->outfp = outfp;
  /* register option -o */
  opto = gt_option_new_string("o", "redirect output to specified file",
                              ofi->output_filename, NULL);
  gt_option_parser_add_option(op, opto);
  /* register option -gzip */
  optgzip = gt_option_new_bool("gzip", "write gzip compressed output file",
                               &ofi->gzip, false);
  gt_option_parser_add_option(op, optgzip);
  /* register option -bzip2 */
  optbzip2 = gt_option_new_bool("bzip2", "write bzip2 compressed output file",
                                &ofi->bzip2, false);
  gt_option_parser_add_option(op, optbzip2);
  /* register option -force */
  optforce = gt_option_new_bool(GT_FORCE_OPT_CSTR,
                                "force writing to output file",
                                &ofi->force, false);
  gt_option_parser_add_option(op, optforce);
  /* options -gzip and -bzip2 exclude each other */
  gt_option_exclude(optgzip, optbzip2);
  /* option implications */
  gt_option_imply(optgzip, opto);
  gt_option_imply(optbzip2, opto);
  gt_option_imply(optforce, opto);
  /* set hook function to determine <outfp> */
  gt_option_parser_register_hook(op, determine_outfp, ofi);
}

void gt_output_file_info_delete(GtOutputFileInfo *ofi)
{
  if (!ofi) return;
  gt_str_delete(ofi->output_filename);
  gt_free(ofi);
}

GtFile* gt_output_file_xopen_forcecheck(const char *path, const char *mode,
                                       bool force, GtError *err)
{
  gt_error_check(err);
  gt_assert(path && mode);
  if (!force && gt_file_exists(path)) {
    gt_error_set(err, "file \"%s\" exists already, use option -%s to overwrite",
                 path, GT_FORCE_OPT_CSTR);
    return NULL;
  }
  return gt_file_xopen(path, mode);
}
