/*
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHAciOEVER RESULTING FROM LOSS OF USE, DATA OR PROFIci, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef ENCSEQ_OPTIONS_H
#define ENCSEQ_OPTIONS_H

#include "core/option_api.h"
#include "core/str_api.h"
#include "core/str_array.h"

typedef struct GtEncseqOptions GtEncseqOptions;

/* Registers the necessary options for creating an encoded sequence in <op>.
   Additionally includes options to specify index name and input files,
   which will be written to <indexname> and <indb>. If not given, index
   names will be generated from the input sequence filename, if there is
   only one.  */
GtEncseqOptions* gt_encseq_options_register_encoding(GtOptionParser *op,
                                                     GtStr *indexname,
                                                     GtStrArray *indb);

/* Adds the `-dir' option to <op>, allowing the selection of a readmode.
   Not all tools support readmodes, so this needs to be added separately.
   The selected value is written to <readmode>. Default is `fwd'. */
void             gt_encseq_options_add_readmode_option(GtOptionParser *op,
                                                       GtStr *readmode);

/* Registers the necessary options for loading an encoded sequence in <op>.
   Returns an info object which can be used to access the options for further
   configuration.
   Additionally includes options to specify an input index name, which will be
   written to <indexname>. */
GtEncseqOptions* gt_encseq_options_register_loading(GtOptionParser *op,
                                                    GtStr *indexname);

void             gt_encseq_options_delete(GtEncseqOptions *oi);

/* XXX: clean this up */
#ifndef GT_ENCSEQ_OPTS_GETTER_DECLS_DEFINED

#define GT_ENCSEQ_OPTS_GETTER_DECL_OPT(VARNAME) \
GtOption* gt_encseq_options_##VARNAME##_option(GtEncseqOptions *i);

#define GT_ENCSEQ_OPTS_GETTER_DECL_VAL(VARNAME, TYPE) \
TYPE gt_encseq_options_##VARNAME##_value(GtEncseqOptions *i);

#define GT_ENCSEQ_OPTS_GETTER_DECL(VARNAME,TYPE) \
GT_ENCSEQ_OPTS_GETTER_DECL_OPT(VARNAME) \
GT_ENCSEQ_OPTS_GETTER_DECL_VAL(VARNAME, TYPE)

#define GT_ENCSEQ_OPTS_GETTER_DECLS_DEFINED
#endif

GT_ENCSEQ_OPTS_GETTER_DECL(db, GtStrArray*);
GT_ENCSEQ_OPTS_GETTER_DECL(des, bool);
GT_ENCSEQ_OPTS_GETTER_DECL(dna, bool);
GT_ENCSEQ_OPTS_GETTER_DECL(lossless, bool);
GT_ENCSEQ_OPTS_GETTER_DECL(md5, bool);
GT_ENCSEQ_OPTS_GETTER_DECL(mirrored, bool);
GT_ENCSEQ_OPTS_GETTER_DECL(plain, bool);
GT_ENCSEQ_OPTS_GETTER_DECL(protein, bool);
GT_ENCSEQ_OPTS_GETTER_DECL(sat, GtStr*);
GT_ENCSEQ_OPTS_GETTER_DECL(sds, bool);
GT_ENCSEQ_OPTS_GETTER_DECL(smap, GtStr*);
GT_ENCSEQ_OPTS_GETTER_DECL(ssp, bool);
GT_ENCSEQ_OPTS_GETTER_DECL(tis, bool);
GT_ENCSEQ_OPTS_GETTER_DECL_OPT(dir);

#endif
