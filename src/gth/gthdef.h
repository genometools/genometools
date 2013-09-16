/*
  Copyright (c) 2003-2007, 2012 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2003-2007       Center for Bioinformatics, University of Hamburg

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

#ifndef GTHDEF_H
#define GTHDEF_H

/* file suffixes for different indices */
#define DNASUFFIX               "dna"
#define POLYASUFFIX             "polya"
#define MAXSUFFIXLEN            14

/* the name of the environment variable containing the path for gth data */
#define GTHDATAENVNAME          "GTHDATADIR"

/* the name of the gth data file path relative to binary */
#define GTHDATADIRNAME          "gthdata"

/* the name of the environment variable used to disable file locking */
#define GTHNOFLOCKENVNAME       "GTHNOFLOCK"

#endif
