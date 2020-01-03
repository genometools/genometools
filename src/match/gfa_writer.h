/*
  Copyright (c) 2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#ifndef GFA_WRITER_H
#define GFA_WRITER_H

#include "core/file.h"
#include "core/encseq.h"
#include "core/error_api.h"

/* The <GtGfaWriter> class allows one to write string graph information in
   the GFA format of SGA (Simpsons and Durbin, 2012). */
typedef struct GtGfaWriter GtGfaWriter;

/* Supported versions of the gfa format. */
typedef enum
{
  GT_GFA_VERSION_1_0,
  GT_GFA_VERSION_2_0
} GtGfaVersion;

/* Creates a new <GtGfaWriter> object, using <file> for output
   and <encseq> as source of information.
   <version> is the GFA version. */
GtGfaWriter* gt_gfa_writer_new(GtFile *file,
                               const GtEncseq *encseq,
                               GtGfaVersion version);

/* Writes the header using the provided information: <minlen> is the
   minimal match length, <inputfilename> is a string to use as input filename,
   <has_containments> shall be true, if containments
   are present, <has_transitives> shall be true, if transitive edges are
   present. Returns 0 on success, -1 on error and sets <err>. */
int           gt_gfa_writer_show_header(GtGfaWriter *aw,
                                        GtUword minlen,
                                        const char *inputfilename,
                                        bool has_containments,
                                        bool has_transitives,
                                        GtError *err);

/* Writes the vertices. Returns 0 on success, -1 on error and sets <err>. */
int           gt_gfa_writer_show_segments(GtGfaWriter *aw,
                                          GtError *err);

/* Writes an edge using Readjoiner SPM information. */
void          gt_spmproc_show_gfa(GtUword suffix_readnum,
                                  GtUword prefix_readnum,
                                  GtUword length,
                                  bool suffixseq_direct,
                                  bool prefixseq_direct,
                                  void *gfa_writer);

/* Deletes a <GtGfa> object. */
void          gt_gfa_writer_delete(GtGfaWriter *aw);

#endif
