/*
  Copyright (c) 2004-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2004-2007 Center for Bioinformatics, University of Hamburg

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

#ifndef GTHXML_H
#define GTHXML_H

#include <stdbool.h>
#include "core/file.h"

/* version for ``regular'' results */
#define GTH_XML_VERSION                    "1.1"

/* version for showing intermediate results */
#define GTH_SPLICED_ALIGNMENT_XML_VERSION  "1.0"

void gth_xml_show_leader(bool intermediate, GtFile *outfp);
void gth_xml_show_trailer(bool intermediate, GtFile *outfp);

#endif
