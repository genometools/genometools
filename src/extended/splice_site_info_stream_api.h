/*
  Copyright (c) 2007-2011 Gordon Gremme <gordon@gremme.org>
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

#ifndef SPLICE_SITE_INFO_STREAM_API_H
#define SPLICE_SITE_INFO_STREAM_API_H

#include <stdio.h>
#include "extended/node_stream_api.h"
#include "extended/region_mapping_api.h"

/* The <GtSpliceSiteInfoStream> is a <GtNodeStream> that gathers splice
   site information from <GtFeatureNode>s. */
typedef struct GtSpliceSiteInfoStream GtSpliceSiteInfoStream;

const GtNodeStreamClass* gt_splice_site_info_stream_class(void);

/* Create a GtSpliceSiteInfoStream, takes ownership of <region_mapping>. */
GtNodeStream* gt_splice_site_info_stream_new(GtNodeStream *in_stream,
                                             GtRegionMapping *region_mapping);
/* Prints splice site information gathered in <ns> to <outfp>. */
bool          gt_splice_site_info_stream_show(GtNodeStream *ns, GtFile *outfp);
/* Returns <true> if an intron has been processed in <ns>, <false> otherwise */
bool          gt_splice_site_info_stream_intron_processed(GtNodeStream *ns);
/* Print information for canonical splice sites as stored in <ns>. */
bool          gt_splice_site_info_stream_show_canonical(GtNodeStream *ns,
                                                        bool show_gc);

#endif
