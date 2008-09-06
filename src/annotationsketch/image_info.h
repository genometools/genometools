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

#ifndef IMAGE_INFO_H
#define IMAGE_INFO_H

#include "annotationsketch/recmap.h"
#include "core/error.h"

typedef struct GT_ImageInfo GT_ImageInfo;

GT_ImageInfo* gt_image_info_new();
void          gt_image_info_delete(GT_ImageInfo*);
unsigned int  gt_image_info_get_height(GT_ImageInfo*);
void          gt_image_info_set_height(GT_ImageInfo*, unsigned int);
/* Takes ownership of <recmap>. */
void          gt_image_info_add_recmap(GT_ImageInfo*, GT_RecMap *recmap);
unsigned long gt_image_info_num_of_recmaps(GT_ImageInfo*);
GT_RecMap*       gt_image_info_get_recmap(GT_ImageInfo*, unsigned long);
void          gt_image_info_get_recmap_ptr(GT_ImageInfo*, GT_RecMap*,
                                           unsigned long);
int           gt_image_info_unit_test(Error*);
#endif
