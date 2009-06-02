/*
  Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#ifndef CUSTOM_TRACK_REP_H
#define CUSTOM_TRACK_REP_H

#include <stdio.h>
#include "annotationsketch/custom_track.h"
#include "core/error.h"
#include "core/range.h"
#include "core/str.h"

typedef int           (*GtCustomTrackRenderFunc)(GtCustomTrack*, GtGraphics*,
                                                 unsigned int, GtRange,
                                                 GtStyle*, GtError*);
typedef unsigned long (*GtCustomTrackGetHeightFunc)(GtCustomTrack*);
typedef const char*   (*GtCustomTrackGetTitleFunc)(GtCustomTrack*);
typedef void          (*GtCustomTrackFreeFunc)(GtCustomTrack*);

typedef struct GtCustomTrackMembers GtCustomTrackMembers;

struct GtCustomTrack {
  const GtCustomTrackClass *c_class;
  GtCustomTrackMembers *pvt;
};

const GtCustomTrackClass* gt_custom_track_class_new(size_t size,
                                          GtCustomTrackRenderFunc render,
                                          GtCustomTrackGetHeightFunc get_height,
                                          GtCustomTrackGetTitleFunc get_title,
                                          GtCustomTrackFreeFunc free);
GtCustomTrack* gt_custom_track_create(const GtCustomTrackClass*);
void*          gt_custom_track_cast(const GtCustomTrackClass*, GtCustomTrack*);

#endif
