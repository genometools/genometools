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

#ifndef ASQG_WRITER_H
#define ASQG_WRITER_H

#include "core/file.h"
#include "core/encseq.h"
#include "core/error.h"

/* The <GtAsqgWriter> class allows one to write string graph information in
   the ASQG format of SGA (Simpsons and Durbin, 2012). */
typedef struct GtAsqgWriter GtAsqgWriter;

/* <GT_ASQG_VERSION> is the supported version of the asqg format. */
#define GT_ASQG_VERSION 1

/* Creates a new <GtAsqgWriter> object, using <file> for output
   and <encseq> as source of information. */
GtAsqgWriter* gt_asqg_writer_new(GtFile *file, const GtEncseq *encseq);

/* Writes the header using the provided information: <minlen> is the
   minimal match length, <inputfilename> is a string to use as input filename,
   <erate> the error rate, <has_containments> shall be true, if containments
   are present, <has_transitives> shall be true, if transitive edges are
   present. Returns 0 on success, -1 on error and sets <err>. */
int           gt_asqg_writer_show_header(GtAsqgWriter *aw,
                                         float erate,
                                         GtUword minlen,
                                         const char *inputfilename,
                                         bool has_containments,
                                         bool has_transitives,
                                         GtError *err);

/* Writes the vertices. Returns 0 on success, -1 on error and sets <err>. */
int           gt_asqg_writer_show_vertices(GtAsqgWriter *aw,
                                           GtError *err);

/* Writes an edge using Readjoiner SPM information. */
void          gt_spmproc_show_asgq(GtUword suffix_readnum,
                                   GtUword prefix_readnum,
                                   GtUword length,
                                   bool suffixseq_direct,
                                   bool prefixseq_direct,
                                   void *asqg_writer);

/* Deletes a <GtAsqg> object. */
void          gt_asqg_writer_delete(GtAsqgWriter *aw);

#endif
