/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtcore/fileutils.h>
#include <libgtcore/outputfile.h>

#define FORCE_OPT_CSTR  "force"

struct OutputFileInfo {
  Str *output_filename;
  bool gzip,
       bzip2,
       force;
  GenFile **outfp;
};

OutputFileInfo* outputfileinfo_new(Env *env)
{
  OutputFileInfo *ofi;
  env_error_check(env);
  ofi = env_ma_malloc(env, sizeof (OutputFileInfo));
  ofi->output_filename = str_new(env);
  return ofi;
}

static int determine_outfp(void *data, Env *env)
{
  OutputFileInfo *ofi = (OutputFileInfo*) data;
  GenFileMode genfilemode;
  int has_err = 0;
  env_error_check(env);
  assert(ofi);
  if (!str_length(ofi->output_filename)) /* no output file given -> use stdin */
    *ofi->outfp = NULL;
  else { /* outputfile given -> create generic file pointer */
    if (!ofi->force && file_exists(str_get(ofi->output_filename))) {
        env_error_set(env, "file \"%s\" exists already, use option -%s to "
                      "overwrite", str_get(ofi->output_filename),
                      FORCE_OPT_CSTR);
        has_err = -1;
    }
    if (!has_err) {
      assert(!(ofi->gzip && ofi->bzip2));
      if (ofi->gzip)
        genfilemode = GFM_GZIP;
      else if (ofi->bzip2)
        genfilemode = GFM_BZIP2;
      else
        genfilemode = GFM_UNCOMPRESSED;
      /* XXX: append .gz/.bz2 to output_filename if necessary */
      *ofi->outfp = genfile_xopen(genfilemode, str_get(ofi->output_filename),
                                  "w", env);
      assert(*ofi->outfp);
    }
  }
  return has_err;
}

void outputfile_register_options(OptionParser *op, GenFile **outfp,
                                 OutputFileInfo *ofi, Env *env)
{
  Option *opto, *optgzip, *optbzip2, *optforce;
  env_error_check(env);
  assert(outfp && ofi);
  ofi->outfp = outfp;
  /* register option -o */
  opto = option_new_string("o", "redirect output to specified file",
                           ofi->output_filename, NULL, env);
  option_parser_add_option(op, opto, env);
  /* register option -gzip */
  optgzip = option_new_bool("gzip", "write gzip compressed output file",
                            &ofi->gzip, false, env);
  option_parser_add_option(op, optgzip, env);
  /* register option -bzip2 */
  optbzip2 = option_new_bool("bzip2", "write bzip2 compressed output file",
                             &ofi->bzip2, false, env);
  option_parser_add_option(op, optbzip2, env);
  /* register option -force */
  optforce = option_new_bool(FORCE_OPT_CSTR, "force writing to output file",
                             &ofi->force, false, env);
  option_parser_add_option(op, optforce, env);
  /* options -gzip and -bzip2 exclude each other */
  option_exclude(optgzip, optbzip2, env);
  /* set hook function to determine <outfp> */
  option_parser_register_hook(op, determine_outfp, ofi, env);
}

void outputfileinfo_delete(OutputFileInfo *ofi, Env *env)
{
  if (!ofi) return;
  str_delete(ofi->output_filename, env);
  env_ma_free(ofi, env);
}
