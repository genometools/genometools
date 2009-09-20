/*
  Copyright (c) 2008-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008      Center for Bioinformatics, University of Hamburg

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

#ifndef XML_INTER_SA_VISITOR_H
#define XML_INTER_SA_VISITOR_H

#include "gth/sa_visitor.h"

/* implements the ``spliced alignment visitor'' interface */
typedef struct GthXMLInterSAVisitor GthXMLInterSAVisitor;

const GthSAVisitorClass* gth_xml_inter_sa_visitor_class(void);
GthSAVisitor*            gth_xml_inter_sa_visitor_new(GthInput*,
                                                      unsigned int indentlevel,
                                                      GtFile *outfp);
void                     gth_xml_inter_sa_visitor_set_outfp(GthSAVisitor*,
                                                            GtFile*);

#endif
