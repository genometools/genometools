/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#ifndef LTRFILEOUT_STREAM_H
#define LTRFILEOUT_STREAM_H

#include "libgtcore/bioseq.h"
#include "libgtext/genome_stream.h"
#include "libgtltr/pbs.h"
#include "libgtltr/ppt.h"
#include "libgtltr/pdom.h"

/* implements the ``genome_stream'' interface */
typedef struct LTRFileOutStream LTRFileOutStream;

const GenomeStreamClass* ltr_fileout_stream_class(void);

GenomeStream* ltr_fileout_stream_new(GenomeStream *in_stream,
                                     Bioseq *bioseq,
                                     FILE *fp,
                                     bool with_metadata,
                                     PPTOptions *ppt_opts,
                                     PBSOptions *pbs_opts,
                                     PdomOptions *pdom_opts,
                                     const char *trnafilename,
                                     const char *seqfilename,
                                     const char *gfffilename);

#endif
