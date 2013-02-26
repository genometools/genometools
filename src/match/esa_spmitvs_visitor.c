/*
  Copyright (c) 2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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
#include "core/class_alloc_lock.h"
#include "core/mathsupport.h"
#include "core/unused_api.h"
#include "esa-seqread.h"
#include "esa-spmitvs.h"
#include "esa-bottomup.h"
#include "esa_visitor_rep.h"
#include "esa_spmitvs_visitor.h"

typedef struct
{
  unsigned long wholeleaf, wholeleafwidth, nowholeleaf, nowholeleafwidth;
} Lcpintervalcount;

struct GtESASpmitvsVisitor {
  const GtESAVisitor parent_instance;
  unsigned long unnecessaryleaves,
                totallength,
                currentleafindex,
                lastwholeleaf,
                maxlen;
  unsigned int prefixlength;
  Lcpintervalcount *wholeleafcount;
  const GtEncseq *encseq;
  GtReadmode readmode;
};

#define gt_esa_spmitvs_visitor_cast(GV)\
        gt_esa_visitor_cast(gt_esa_spmitvs_visitor_class(), GV)

static bool gt_esa_spmitvs_visitor_iswholeleaf(const GtEncseq *encseq,
                                               GtReadmode readmode,
                                               unsigned long leafnumber)
{
  return (leafnumber > 0)
    ? gt_encseq_position_is_separator(encseq,leafnumber - 1,readmode)
    : true;
}

static int gt_esa_spmitvs_visitor_processleafedge(GtESAVisitor *ev,
                                                  GT_UNUSED bool firstsucc,
                                                  unsigned long fd,
                                                  GT_UNUSED unsigned long flb,
                                                  GT_UNUSED
                                                         GtESAVisitorInfo *info,
                                                  unsigned long leafnumber,
                                                  GT_UNUSED GtError *err)

{
  GtESASpmitvsVisitor *esv;
  gt_assert(ev);
  esv = gt_esa_spmitvs_visitor_cast(ev);

  if (gt_esa_spmitvs_visitor_iswholeleaf(esv->encseq,esv->readmode,leafnumber))
  {
    gt_assert(esv->currentleafindex != esv->totallength);
    esv->lastwholeleaf = esv->currentleafindex;
  } else
  {
    if (leafnumber + fd < esv->totallength &&
        !gt_encseq_position_is_separator(esv->encseq,
                                         leafnumber + fd,
                                         esv->readmode))
    {
      esv->unnecessaryleaves++;
    }
  }
  esv->currentleafindex++;
  return 0;
}

static int gt_esa_spmitvs_visitor_processbranchingedge(GtESAVisitor *ev,
                                                    GT_UNUSED bool firstsucc,
                                                    unsigned long fd,
                                                    GT_UNUSED unsigned long flb,
                                                    GT_UNUSED
                                                        GtESAVisitorInfo *finfo,
                                                    unsigned long sd,
                                                    unsigned long slb,
                                                    unsigned long srb,
                                                    GT_UNUSED
                                                        GtESAVisitorInfo *sinfo,
                                                    GT_UNUSED GtError *err)
{
  GtESASpmitvsVisitor *esv;
  unsigned long idx;
  gt_assert(ev);
  esv = gt_esa_spmitvs_visitor_cast(ev);

  for (idx=fd+1; idx<sd; idx++)
  {
    gt_assert(idx <= esv->maxlen);
    if (esv->lastwholeleaf != esv->totallength &&
        esv->lastwholeleaf >= slb)
    {
      gt_assert(esv->lastwholeleaf <= srb);
      esv->wholeleafcount[idx].wholeleaf++;
      esv->wholeleafcount[idx].wholeleafwidth += (srb - slb + 1);
    } else
    {
      esv->wholeleafcount[idx].nowholeleaf++;
      esv->wholeleafcount[idx].nowholeleafwidth += (srb - slb + 1);
    }
  }
  return 0;
}

static int gt_esa_spmitvs_visitor_processlcpinterval(GtESAVisitor *ev,
                                                     unsigned long lcp,
                                                     unsigned long lb,
                                                     unsigned long rb,
                                                     GT_UNUSED
                                                         GtESAVisitorInfo *info,
                                                     GT_UNUSED GtError *err)
{
  GtESASpmitvsVisitor *esv;
  gt_assert(ev);
  esv = gt_esa_spmitvs_visitor_cast(ev);

  if (esv->lastwholeleaf != esv->totallength &&
      esv->lastwholeleaf >= lb)
  {
    gt_assert(lcp <= (unsigned long) esv->maxlen);
    gt_assert(esv->lastwholeleaf <= rb);
    esv->wholeleafcount[lcp].wholeleaf++;
    esv->wholeleafcount[lcp].wholeleafwidth += (rb - lb + 1);
  } else
  {
    esv->wholeleafcount[lcp].nowholeleaf++;
    esv->wholeleafcount[lcp].nowholeleafwidth += (rb - lb + 1);
  }
  return 0;
}

void gt_esa_spmitvs_visitor_delete(GtESAVisitor *ev)
{
  GtESASpmitvsVisitor *esv;
  gt_assert(ev);
  esv = gt_esa_spmitvs_visitor_cast(ev);
  gt_free(esv->wholeleafcount);
}

const GtESAVisitorClass* gt_esa_spmitvs_visitor_class()
{
  static const GtESAVisitorClass *esc = NULL;
  gt_class_alloc_lock_enter();
  if (!esc) {
    esc = gt_esa_visitor_class_new(sizeof (GtESASpmitvsVisitor),
                                   gt_esa_spmitvs_visitor_delete,
                                   gt_esa_spmitvs_visitor_processleafedge,
                                   gt_esa_spmitvs_visitor_processbranchingedge,
                                   gt_esa_spmitvs_visitor_processlcpinterval,
                                   NULL,
                                   NULL);
  }
  gt_class_alloc_lock_leave();
  return esc;
}

GtESAVisitor* gt_esa_spmitvs_visitor_new(const GtEncseq *encseq,
                                         GtReadmode readmode,
                                         unsigned int prefixlength,
                                         GT_UNUSED GtError *err)
{
  GtESAVisitor *ev = gt_esa_visitor_create(gt_esa_spmitvs_visitor_class());
  GtESASpmitvsVisitor *esv = gt_esa_spmitvs_visitor_cast(ev);

  esv->encseq = encseq;
  esv->readmode = readmode;
  esv->maxlen = gt_encseq_max_seq_length(encseq);
  esv->unnecessaryleaves = 0;
  esv->totallength = gt_encseq_total_length(encseq);
  esv->currentleafindex = 0;
  esv->lastwholeleaf = esv->totallength; /* undefined */
  esv->prefixlength = prefixlength;
  esv->wholeleafcount = gt_malloc(sizeof (*esv->wholeleafcount) *
                                 (esv->maxlen+1));
  memset(esv->wholeleafcount, 0,
         sizeof (*esv->wholeleafcount) * (esv->maxlen+1));

  return ev;
}

void gt_esa_spmitvs_visitor_print_results(GtESASpmitvsVisitor *esv,
                                          unsigned long nonspecials)
{
  unsigned long idx;
  printf("unnecessaryleaves=%lu (%.2f)\n",
         esv->unnecessaryleaves,
         (double) esv->unnecessaryleaves/nonspecials);
  for (idx = 0; idx<= esv->maxlen; idx++)
  {
    if (esv->wholeleafcount[idx].wholeleaf != 0 ||
        esv->wholeleafcount[idx].nowholeleaf != 0)
    {
      printf("wholeleaf[%lu]:num=%lu (%.2f), ",idx,
             esv->wholeleafcount[idx].wholeleaf,
             (double) esv->wholeleafcount[idx].wholeleaf/
                      (esv->wholeleafcount[idx].wholeleaf+
                       esv->wholeleafcount[idx].nowholeleaf));
      printf("width=%lu (%.2f)\n",
              esv->wholeleafcount[idx].wholeleafwidth,
              (double) esv->wholeleafcount[idx].wholeleafwidth/
                       esv->totallength);
    }
  }
}
