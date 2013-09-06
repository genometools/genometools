/*
  Copyright (c) 2008-2009 Gordon Gremme <gordon@gremme.org>
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

#ifndef XML_PGL_VISITOR_H
#define XML_PGL_VISITOR_H

#include "gth/input.h"
#include "gth/pgl_visitor.h"

/* implements the ``spliced alignment visitor'' interface */
typedef struct GthXMLPGLVisitor GthXMLPGLVisitor;

const GthPGLVisitorClass* gth_xml_pgl_visitor_class(void);
GthPGLVisitor*            gth_xml_pgl_visitor_new(GthInput*,
                                                  GtUword
                                                  translationtable,
                                                  unsigned int indentlevel,
                                                  GthOutput*);

#endif
