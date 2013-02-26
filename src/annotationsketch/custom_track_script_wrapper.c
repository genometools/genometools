  /*
  Copyright (c) 2008-2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2009 Center for Bioinformatics, University of Hamburg

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

#include "annotationsketch/custom_track_script_wrapper.h"
#include "annotationsketch/custom_track_rep.h"
#include "core/class_alloc_lock.h"

struct GtCustomTrackScriptWrapper {
  const GtCustomTrack parent_instance;
  GtCtScriptRenderFunc render_func;
  GtCtScriptGetHeightFunc get_height_func;
  GtCtScriptGetTitleFunc get_title_func;
  GtCtScriptFreeFunc free_func;
  GtStr *title_buffer;
};

#define gt_custom_track_script_wrapper_cast(ct) \
        gt_custom_track_cast(gt_custom_track_script_wrapper_class(), ct)

int gt_custom_track_script_wrapper_sketch(GtCustomTrack *ct,
                                          GtGraphics *graphics,
                                          unsigned int start_ypos,
                                          GtRange viewrange,
                                          GtStyle *style,
                                          GtError *err)
{
  GtCustomTrackScriptWrapper *cte;
  cte = gt_custom_track_script_wrapper_cast(ct);
  return cte->render_func(graphics, start_ypos, &viewrange, style, err);
}

unsigned long gt_custom_track_script_wrapper_get_height(GtCustomTrack *ct)
{
  GtCustomTrackScriptWrapper *cte;
  cte = gt_custom_track_script_wrapper_cast(ct);
  return cte->get_height_func(NULL);
}

const char* gt_custom_track_script_wrapper_get_title(GtCustomTrack *ct)
{
  GtCustomTrackScriptWrapper *cte;
  cte = gt_custom_track_script_wrapper_cast(ct);
  gt_str_reset(cte->title_buffer);
  cte->get_title_func(NULL, cte->title_buffer);
  return gt_str_get(cte->title_buffer);
}

void gt_custom_track_script_wrapper_delete(GtCustomTrack *ct)
{
  GtCustomTrackScriptWrapper *cte;
  if (!ct) return;
  cte = gt_custom_track_script_wrapper_cast(ct);
  cte->free_func(NULL);
  gt_str_delete(cte->title_buffer);
}

const GtCustomTrackClass* gt_custom_track_script_wrapper_class(void)
{
  static const GtCustomTrackClass *ctc = NULL;
  gt_class_alloc_lock_enter();
  if (!ctc)
  {
    ctc = gt_custom_track_class_new(sizeof (GtCustomTrackScriptWrapper),
                                    gt_custom_track_script_wrapper_sketch,
                                    gt_custom_track_script_wrapper_get_height,
                                    gt_custom_track_script_wrapper_get_title,
                                    gt_custom_track_script_wrapper_delete);
  }
  gt_class_alloc_lock_leave();
  return ctc;
}

GtCustomTrack* gt_custom_track_script_wrapper_new(GtCtScriptRenderFunc
                                                             render_func,
                                                  GtCtScriptGetHeightFunc
                                                             get_height_func,
                                                  GtCtScriptGetTitleFunc
                                                             get_title_func,
                                                  GtCtScriptFreeFunc
                                                             free_func)
{
  GtCustomTrackScriptWrapper *cte;
  GtCustomTrack *ct;
  ct = gt_custom_track_create(gt_custom_track_script_wrapper_class());
  cte = gt_custom_track_script_wrapper_cast(ct);
  cte->render_func = render_func;
  cte->get_height_func = get_height_func;
  cte->get_title_func = get_title_func;
  cte->free_func = free_func;
  cte->title_buffer = gt_str_new();
  return ct;
}
