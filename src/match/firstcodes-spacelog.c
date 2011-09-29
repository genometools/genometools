#include "core/ma.h"
#include "core/log_api.h"
#include "core/spacecalc.h"
#include "firstcodes-spacelog.h"

struct GtFirstcodesspacelog
{
  size_t workspace, splitspace;
};

GtFirstcodesspacelog *gt_firstcodes_spacelog_new(void)
{
  GtFirstcodesspacelog *fcsl;

  fcsl = gt_malloc(sizeof (GtFirstcodesspacelog));
  fcsl->workspace = fcsl->splitspace = 0;
  return fcsl;
}

void gt_firstcodes_spacelog_delete(GtFirstcodesspacelog *fcsl)
{
  gt_free(fcsl);
}

size_t gt_firstcodes_spacelog_total(GtFirstcodesspacelog *fcsl)
{
  return fcsl->workspace + fcsl->splitspace;
}

void gt_firstcodes_spacelog_add(GtFirstcodesspacelog *fcsl,
                                int line,
                                const char *filename,
                                bool add,
                                const char *kind,
                                bool addtowork,
                                size_t workspace)
{
  if (addtowork)
  {
    if (add)
    {
      fcsl->workspace += workspace;
    } else
    {
      gt_assert(fcsl->workspace >= workspace);
      fcsl->workspace -= workspace;
    }
  } else
  {
    if (add)
    {
      fcsl->splitspace += workspace;
    } else
    {
      gt_assert(fcsl->workspace >= workspace);
      fcsl->splitspace -= workspace;
    }
  }
  gt_log_log("file %s, line %d: %s %.2f MB for %s to %space; work=%.2f, "
             "split=%.2f,all=%.2f MB",
             filename,
             line,
             add ? "add" : "delete",
             GT_MEGABYTES(workspace),
             kind,
             addtowork ? "work" : "split",
             GT_MEGABYTES(fcsl->workspace),
             GT_MEGABYTES(fcsl->splitspace),
             GT_MEGABYTES(fcsl->workspace+fcsl->splitspace));
  /*
  gt_logger_log(logger,"current space peak %.2f",
                GT_MEGABYTES(gt_spacepeak_get_space_peak()));
  */
}
