/*
  Copyright (c) 2007-2012 Gordon Gremme <gordon@gremme.org>
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

#ifndef SEQID2FILE_API_H
#define SEQID2FILE_API_H

#include "core/option_api.h"
#include "extended/region_mapping_api.h"

/* The <GtSeqid2FileInfo> class represents the state of a sequence source
   statement as given by the -seqfile, -seqfiles, -matchdesc, -usedesc and
   -regionmapping options. */
typedef struct GtSeqid2FileInfo GtSeqid2FileInfo;

/* Create a new <GtSeqid2FileInfo> object. */
GtSeqid2FileInfo* gt_seqid2file_info_new(void);
/* Create a new <GtSeqid2FileInfo> object. */
void              gt_seqid2file_info_delete(GtSeqid2FileInfo *);

/* SeqID2File module */

/* Add the options -seqfile, -seqfiles, -matchdesc, -usedesc and
   -regionmapping to the given <option_parser>. */
void              gt_seqid2file_register_options(GtOptionParser *option_parser,
                                                 GtSeqid2FileInfo *s2fi);
/* Add the options -seqfile, -seqfiles, -matchdesc, -usedesc and
   -regionmapping to the given <option_parser>. If <mandatory> is set,
   either option -seqfile, -seqfiles or -regionmapping is mandatory.
   If <debug> is set, then the options are marked as development options. */
void              gt_seqid2file_register_options_ext(GtOptionParser
                                                     *option_parser,
                                                     GtSeqid2FileInfo *s2fi,
                                                     bool mandatory,
                                                     bool debug);
/* Returns TRUE if any of the options -seqfile, -seqfiles, -matchdesc,
   -usedesc or -regionmapping stored in <s2fi> has been specified and
   given a parameter. */
bool              gt_seqid2file_option_used(GtSeqid2FileInfo *s2fi);
/* Returns a <GtRegionMapping> based on the <s2fi>. <NULL> will be returned
   on error, and <err> will be set accordingly. */
GtRegionMapping*  gt_seqid2file_region_mapping_new(GtSeqid2FileInfo *s2fi,
                                                   GtError *err);

#endif
