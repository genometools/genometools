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

#include <math.h>
#include "core/unused_api.h"
#include "core/readmode_api.h"
#include "core/error_api.h"
#include "esa-bottomup.h"
#include "esa-spmsk.h"

typedef struct
{
  unsigned long firstinW;
} GtSpmsk_info;

struct GtSpmsk_state /* global information */
{
  unsigned long totallength,
                minmatchlength;
  const GtEncseq *encseq;
  GtReadmode readmode;
  GtSpmsk_info *info;
};

static GtBUinfo *gt_spmsk_allocatestackinfo(GtBUstate *bustate)
{
  return (GtBUinfo *) (((GtSpmsk_state *) bustate)->info);
}

static int spmsk_processleafedge(GT_UNUSED bool firstsucc,
                                 GT_UNUSED unsigned long fd,
                                 GT_UNUSED unsigned long flb,
                                 GT_UNUSED GtBUinfo *info,
                                 GT_UNUSED unsigned long leafnumber,
                                 GtBUstate *bustate,
                                 GT_UNUSED GtError *err)

{
  GT_UNUSED GtSpmsk_state *spmsk_state = (GtSpmsk_state *) bustate;

  return 0;
}

static int spmsk_processbranchingedge(GT_UNUSED bool firstsucc,
                                      GT_UNUSED unsigned long fd,
                                      GT_UNUSED unsigned long flb,
                                      GT_UNUSED GtBUinfo *finfo,
                                      GT_UNUSED unsigned long sd,
                                      GT_UNUSED unsigned long slb,
                                      GT_UNUSED unsigned long srb,
                                      GT_UNUSED GtBUinfo *sinfo,
                                      GtBUstate *bustate,
                                      GT_UNUSED GtError *err)
{
  GT_UNUSED GtSpmsk_state *spmsk_state = (GtSpmsk_state *) bustate;

  return 0;
}

static int spmsk_processlcpinterval(GT_UNUSED unsigned long lcp,
                                    GT_UNUSED unsigned long lb,
                                    GT_UNUSED unsigned long rb,
                                    GT_UNUSED GtBUinfo *info,
                                    GtBUstate *bustate,
                                    GT_UNUSED GtError *err)
{
  GT_UNUSED GtSpmsk_state *spmsk_state = (GtSpmsk_state *) bustate;

  return 0;
}

GtSpmsk_state *gt_spmsk_new(const GtEncseq *encseq,
                            GtReadmode readmode,
                            unsigned long minmatchlength)
{
  GtSpmsk_state *state = gt_malloc(sizeof (*state));

  state->encseq = encseq;
  state->readmode = readmode;
  state->totallength = gt_encseq_total_length(encseq);
  state->minmatchlength = minmatchlength;
  return state;
}

void gt_spmsk_delete(GtSpmsk_state *state)
{
  gt_free(state);
}

int gt_spmsk_process(GtSpmsk_state *state,
                        const unsigned long *suftab_bucket,
                        const uint16_t *lcptab_bucket,
                        unsigned long nonspecials,
                        GtError *err)
{
  if (gt_esa_bottomup_RAM(suftab_bucket,
                          lcptab_bucket,
                          nonspecials,
                          gt_spmsk_allocatestackinfo,
                          NULL,
                          spmsk_processleafedge,
                          spmsk_processbranchingedge,
                          spmsk_processlcpinterval,
                          (GtBUstate *) &state,
                          err) != 0)
  {
    return -1;
  }
  return 0;
}
