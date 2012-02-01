/*
  Copyright (c) 2010-2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2010-2011 Center for Bioinformatics, University of Hamburg

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

#ifndef RDJ_CNTLIST_H
#define RDJ_CNTLIST_H

/*
 * Two formats for contained reads lists are defined here.
 * One is a text-based one, while the binary one is a copy
 * of the bitsequence preceded by a short header.
 */

#include "core/intbits.h"           /* GtBitsequence */
#include "core/error_api.h"         /* GtError */

/* writes cntlist to file */
int gt_cntlist_show(GtBitsequence *cntlist, unsigned long nofreads,
    const char *path, bool binary, GtError *err);

/* prepare a file for output of cntlist in bin format */
void gt_cntlist_write_bin_header(unsigned long nofreads, FILE *file);

/* parses a cntlist file, format is automatically recognized;
 * if alloc_cntlist is true, **cntlist is allocated, otherwise
 * **cntlist must be a valid bit sequence with the correct size */
int gt_cntlist_parse(const char *filename, bool alloc_cntlist,
    GtBitsequence **cntlist, unsigned long *nofreads, GtError *err);

/* variant of parse: parses, checks that the number of reads is the one
 * specified, returns number of cnt reads, exits on error */
unsigned long gt_cntlist_xload(const char *filename, GtBitsequence **cntlist,
    unsigned long expected_nofreads);

/* counts the number of bit set in bitsequence cntlist */
unsigned long gt_cntlist_count(const GtBitsequence *cntlist,
    unsigned long nofreads);

#endif
