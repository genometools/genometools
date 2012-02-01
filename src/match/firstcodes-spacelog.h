/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#ifndef FIRSTCODES_SPACELOG_H
#define FIRSTCODES_SPACELOG_H

#include <stdio.h>
#include <stdbool.h>

typedef struct GtFirstcodesspacelog GtFirstcodesspacelog;

GtFirstcodesspacelog *gt_firstcodes_spacelog_new(void);

void gt_firstcodes_spacelog_delete(GtFirstcodesspacelog *fcsl);

size_t gt_firstcodes_spacelog_total(const GtFirstcodesspacelog *fcsl);

size_t gt_firstcodes_spacelog_workspace(const GtFirstcodesspacelog *fcsl);

void gt_firstcodes_spacelog_start_diff(GtFirstcodesspacelog *fcsl);

void gt_firstcodes_spacelog_stop_diff(GtFirstcodesspacelog *fcsl);

size_t gt_firstcodes_spacelog_peak(const GtFirstcodesspacelog *fcsl);

bool gt_firstcodes_spacelog_showentries(FILE *fp,
                                        const GtFirstcodesspacelog *fcsl);

void gt_firstcodes_spacelog_add(GtFirstcodesspacelog *fcsl,
                                int line,
                                const char *filename,
                                bool add,
                                const char *title,
                                bool work,
                                size_t size);

#define GT_FCI_ADDWORKSPACE(FCSL,TAB,SIZE)\
        gt_firstcodes_spacelog_add(FCSL,__LINE__,__FILE__,true,TAB,true,SIZE)

#define GT_FCI_ADDSPLITSPACE(FCSL,TAB,SIZE)\
        gt_firstcodes_spacelog_add(FCSL,__LINE__,__FILE__,true,TAB,false,SIZE)

#define GT_FCI_SUBTRACTWORKSPACE(FCSL,TAB)\
        gt_firstcodes_spacelog_add(FCSL,__LINE__,__FILE__,false,TAB,true,0)

#define GT_FCI_SUBTRACTSPLITSPACE(FCSL,TAB)\
        gt_firstcodes_spacelog_add(FCSL,__LINE__,__FILE__,false,TAB,false,0)

#define GT_FCI_SUBTRACTADDWORKSPACE(FCSL,TAB,SIZE)\
        gt_firstcodes_spacelog_add(FCSL,__LINE__,__FILE__,false,TAB,true,SIZE)

#define GT_FCI_SUBTRACTADDSPLITSPACE(FCSL,TAB,SIZE)\
        gt_firstcodes_spacelog_add(FCSL,__LINE__,__FILE__,false,TAB,false,SIZE)

#endif
