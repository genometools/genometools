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

#ifndef INDEX_OPTIONS_H
#define INDEX_OPTIONS_H

#include "core/encseq_options.h"
#include "core/option_api.h"
#include "core/readmode.h"
#include "match/sfx-strategy.h"

typedef struct GtIndexOptions GtIndexOptions;

#define GT_PREFIXLENGTH_AUTOMATIC 0
#define GT_MAXDEPTH_AUTOMATIC     0

/* This module encapsulates the registration of options for index generation. */

GtIndexOptions* gt_index_options_register_esa(GtOptionParser *op,
                                              GtEncseqOptions *encopts);
/* this will only options concerning the building of the index, used for
   genomediff, for example, where no output is produced. */
GtIndexOptions* gt_index_options_register_esa_noout(GtOptionParser *op);
GtIndexOptions* gt_index_options_register_packedidx(GtOptionParser *op,
                                                    GtStr *indexname,
                                                    GtEncseqOptions *encopts);

void            gt_index_options_delete(GtIndexOptions *oi);

int             gt_parse_algbounds(Sfxstrategy *sfxstrategy,
                                   const GtStrArray *algbounds,
                                   GtError *err);

/* XXX: clean this up */
#ifndef GT_INDEX_OPTS_GETTER_DECLS_DEFINED
#define GT_INDEX_OPTS_GETTER_DECL_OPT(VARNAME) \
GtOption* gt_index_options_##VARNAME##_option(GtIndexOptions *i);
#define GT_INDEX_OPTS_GETTER_DECL_VAL(VARNAME, TYPE) \
TYPE gt_index_options_##VARNAME##_value(GtIndexOptions *i);
#define GT_INDEX_OPTS_GETTER_DECL(VARNAME,TYPE) \
GT_INDEX_OPTS_GETTER_DECL_OPT(VARNAME) \
GT_INDEX_OPTS_GETTER_DECL_VAL(VARNAME, TYPE)
#define GT_INDEX_OPTS_GETTER_DECLS_DEFINED
#endif

GT_INDEX_OPTS_GETTER_DECL(algbounds, GtStrArray*);
GT_INDEX_OPTS_GETTER_DECL(outbcktab, bool);
GT_INDEX_OPTS_GETTER_DECL(outbwttab, bool);
GT_INDEX_OPTS_GETTER_DECL(outkystab, bool);
GT_INDEX_OPTS_GETTER_DECL(outlcptab, bool);
GT_INDEX_OPTS_GETTER_DECL(outsuftab, bool);
GT_INDEX_OPTS_GETTER_DECL(prefixlength, unsigned int);
GT_INDEX_OPTS_GETTER_DECL_OPT(spmopt);
GT_INDEX_OPTS_GETTER_DECL_VAL(bwtIdxParams, struct bwtOptions);
GT_INDEX_OPTS_GETTER_DECL_VAL(lcpdist, bool);
GT_INDEX_OPTS_GETTER_DECL_VAL(maximumspace, unsigned long);
GT_INDEX_OPTS_GETTER_DECL_VAL(numofparts, unsigned int);
GT_INDEX_OPTS_GETTER_DECL_VAL(outkyssort, bool);
GT_INDEX_OPTS_GETTER_DECL_VAL(readmode, GtReadmode);
GT_INDEX_OPTS_GETTER_DECL_VAL(sfxstrategy, Sfxstrategy);
GT_INDEX_OPTS_GETTER_DECL_VAL(swallow_tail, bool);

#endif
