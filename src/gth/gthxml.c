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

#include "gth/gthxml.h"

void gth_xml_show_leader(bool intermediate, GtFile *outfp)
{
  gt_file_xprintf(outfp,
                     "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
  if (intermediate) {
    gt_file_xprintf(outfp, "<SplicedAlignment xmlns="
                    "\"http://www.GenomeThreader.org/SplicedAlignment/\" "
                    "GTH_spliced_alignment_XML_version=\"%s\">\n",
                    GTH_SPLICED_ALIGNMENT_XML_VERSION);
  }
  else {
    gt_file_xprintf(outfp, "<GTH_output xmlns="
                     "\"http://www.genomethreader.org/GTH_output/\" "
                     "GTH_XML_version=\"%s\">\n", GTH_XML_VERSION);
  }
}

void gth_xml_show_trailer(bool intermediate, GtFile *outfp)
{
  if (intermediate)
    gt_file_xprintf(outfp, "</SplicedAlignment>\n");
  else
    gt_file_xprintf(outfp, "</GTH_output>\n");
}
