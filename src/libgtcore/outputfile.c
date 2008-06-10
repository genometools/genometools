/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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
#include "libgtcore/fileutils.h"
#include "libgtcore/ma.h"
#include "libgtcore/outputfile.h"
#include "libgtcore/warning.h"

struct OutputFileInfo {
  Str *output_filename;
  bool gzip,
       bzip2,
       force;
  GenFile **outfp;
};

OutputFileInfo* outputfileinfo_new(void)
{
  OutputFileInfo *ofi;
  ofi = ma_malloc(sizeof (OutputFileInfo));
  ofi->output_filename = str_new();
  return ofi;
}

static int determine_outfp(void *data, Error *err)
{
  OutputFileInfo *ofi = (OutputFileInfo*) data;
  GenFileMode genfilemode;
  int had_err = 0;
  error_check(err);
  assert(ofi);
  if (!str_length(ofi->output_filename)) /* no output file given -> use stdin */
    *ofi->outfp = NULL;
  else { /* outputfile given -> create generic file pointer */
    assert(!(ofi->gzip && ofi->bzip2));
    if (ofi->gzip)
      genfilemode = GFM_GZIP;
    else if (ofi->bzip2)
      genfilemode = GFM_BZIP2;
    else
      genfilemode = GFM_UNCOMPRESSED;
    if (genfilemode != GFM_UNCOMPRESSED &&
        strcmp(str_get(ofi->output_filename) +
               str_length(ofi->output_filename) -
               strlen(genfilemode_suffix(genfilemode)),
               genfilemode_suffix(genfilemode))) {
      warning("output file '%s' doesn't have correct suffix '%s', appending "
              "it", str_get(ofi->output_filename),
              genfilemode_suffix(genfilemode));
      str_append_cstr(ofi->output_filename, genfilemode_suffix(genfilemode));
    }
    if (!ofi->force && file_exists(str_get(ofi->output_filename))) {
        error_set(err, "file \"%s\" exists already, use option -%s to "
                  "overwrite", str_get(ofi->output_filename), FORCE_OPT_CSTR);
        had_err = -1;
    }
    if (!had_err) {
      *ofi->outfp = genfile_xopen_w_gfmode(genfilemode,
                                           str_get(ofi->output_filename), "w");
      assert(*ofi->outfp);
    }
  }
  return had_err;
}

void outputfile_register_options(OptionParser *op, GenFile **outfp,
                                 OutputFileInfo *ofi)
{
  Option *opto, *optgzip, *optbzip2, *optforce;
  assert(outfp && ofi);
  ofi->outfp = outfp;
  /* register option -o */
  opto = option_new_string("o", "redirect output to specified file",
                           ofi->output_filename, NULL);
  option_parser_add_option(op, opto);
  /* register option -gzip */
  optgzip = option_new_bool("gzip", "write gzip compressed output file",
                            &ofi->gzip, false);
  option_parser_add_option(op, optgzip);
  /* register option -bzip2 */
  optbzip2 = option_new_bool("bzip2", "write bzip2 compressed output file",
                             &ofi->bzip2, false);
  option_parser_add_option(op, optbzip2);
  /* register option -force */
  optforce = option_new_bool(FORCE_OPT_CSTR, "force writing to output file",
                             &ofi->force, false);
  option_parser_add_option(op, optforce);
  /* options -gzip and -bzip2 exclude each other */
  option_exclude(optgzip, optbzip2);
  /* option implications */
  option_imply(optgzip, opto);
  option_imply(optbzip2, opto);
  option_imply(optforce, opto);
  /* set hook function to determine <outfp> */
  option_parser_register_hook(op, determine_outfp, ofi);
}

void outputfileinfo_delete(OutputFileInfo *ofi)
{
  if (!ofi) return;
  str_delete(ofi->output_filename);
  ma_free(ofi);
}
