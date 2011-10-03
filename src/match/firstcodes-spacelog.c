#include <string.h>
#include "core/ma.h"
#include "core/fa.h"
#include "core/log_api.h"
#include "core/spacecalc.h"
#include "core/spacepeak.h"
#include "core/mathsupport.h"
#include "firstcodes-spacelog.h"

typedef struct
{
  const char *filename;
  int line;
  const char *title;
  size_t size;
  bool work;
} GtFirstcodespacelogentry;

struct GtFirstcodesspacelog
{
  size_t workspace,
         splitspace,
         spacepeak;
  double max_percent_difference;
  bool calc_difference;
  GtFirstcodespacelogentry *entries;
  unsigned long nextfree, allocated;
};

GtFirstcodesspacelog *gt_firstcodes_spacelog_new(void)
{
  GtFirstcodesspacelog *fcsl;

  fcsl = gt_malloc(sizeof (GtFirstcodesspacelog));
  fcsl->workspace = fcsl->splitspace = 0;
  fcsl->nextfree = fcsl->allocated = 0;
  fcsl->spacepeak = 0;
  fcsl->max_percent_difference = 0.0;
  fcsl->calc_difference = false;
  fcsl->entries = NULL;
  return fcsl;
}

void gt_firstcodes_spacelog_start_diff(GtFirstcodesspacelog *fcsl)
{
  fcsl->calc_difference = (gt_ma_enabled() && gt_fa_enabled()) ? true : false;
}

void gt_firstcodes_spacelog_stop_diff(GtFirstcodesspacelog *fcsl)
{
  fcsl->calc_difference = false;
}

size_t gt_firstcodes_spacelog_total(const GtFirstcodesspacelog *fcsl)
{
  return fcsl->workspace + fcsl->splitspace;
}

size_t gt_firstcodes_spacelog_workspace(const GtFirstcodesspacelog *fcsl)
{
  return fcsl->workspace;
}

size_t gt_firstcodes_spacelog_peak(const GtFirstcodesspacelog *fcsl)
{
  return fcsl->spacepeak;
}

static GtFirstcodespacelogentry *gt_spacelog_find(GtFirstcodesspacelog *fcsl,
                                                  const char *title)
{
  unsigned long idx;

  for (idx = 0; idx < fcsl->nextfree; idx++)
  {
    if (strcmp(fcsl->entries[idx].title,title) == 0)
    {
      return fcsl->entries + idx;
    }
  }
  return NULL;
}

static void gt_spacelog_updateaddentry(GtFirstcodespacelogentry *entry,
                                       const char *filename,int line,
                                       const char *title,size_t size,
                                       bool work)
{
  entry->filename = filename;
  entry->line = line;
  entry->title = title;
  entry->size = size;
  entry->work = work;
}

static void gt_spacelog_addentry(GtFirstcodesspacelog *fcsl,
                                 const char *filename,int line,
                                 const char *title,size_t size,
                                 bool work)
{
  gt_assert(fcsl->nextfree <= fcsl->allocated);
  if (fcsl->nextfree == fcsl->allocated)
  {
    fcsl->allocated += 16UL;
    fcsl->entries = gt_realloc(fcsl->entries,sizeof (*fcsl->entries) *
                                             fcsl->allocated);
  }
  gt_spacelog_updateaddentry(fcsl->entries + fcsl->nextfree,
                             filename,line,title,size,work);
  fcsl->nextfree++;
}

bool gt_firstcodes_spacelog_showentries(FILE *fp,
                                        const GtFirstcodesspacelog *fcsl)
{
  unsigned long idx;
  bool foundnonempty = false;

  for (idx = 0; idx < fcsl->nextfree; idx++)
  {
    if (fcsl->entries[idx].size > 0)
    {
      fprintf(fp,"%s %d %s %s %lu\n",
               fcsl->entries[idx].filename,
               fcsl->entries[idx].line,
               fcsl->entries[idx].title,
               fcsl->entries[idx].work ? "work" : "split",
               (unsigned long) fcsl->entries[idx].size);
      foundnonempty = true;
    }
  }
  return foundnonempty;
}

void gt_firstcodes_spacelog_delete(GtFirstcodesspacelog *fcsl)
{
  if (fcsl != NULL)
  {
#ifdef SKDEBUG
    if (gt_firstcodes_spacelog_showentries(stderr,fcsl))
    {
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
#endif
    gt_log_log("maximal difference between estimated and real space = %.2f%%",
               fcsl->max_percent_difference);
    gt_free(fcsl->entries);
    gt_free(fcsl);
  }
}

static void gt_firstcodes_subtract_error(const char *title,
                                         const char *filename,
                                         int line,
                                         size_t size,
                                         bool work,
                                         size_t sumspace)
{
  fprintf(stderr,"for title \"%s\" (from file %s, line %d) "
                 "in spacelog entries: "
                 "size=%lu > %lu=%sspace\n",title,filename,line,
                 (unsigned long) size,(unsigned long) sumspace,
                 work ? "work" : "split");
}

static void gt_firstcodes_updatemax(GtFirstcodesspacelog *fcsl)
{
  if (fcsl->workspace + fcsl->splitspace > fcsl->spacepeak)
  {
    fcsl->spacepeak = fcsl->workspace + fcsl->splitspace;
    gt_log_log("update spacepeak to %.2f (%lu)",
               GT_MEGABYTES(fcsl->spacepeak),
               (unsigned long) fcsl->spacepeak);
  }
}

void gt_firstcodes_spacelog_add(GtFirstcodesspacelog *fcsl,
                                int line,
                                const char *filename,
                                bool add,
                                const char *title,
                                bool work,
                                size_t size)
{
  GtFirstcodespacelogentry *entry;
  size_t logspace;

  if (add)
  {
    entry = gt_spacelog_find(fcsl,title);
    if (entry != NULL)
    {
      if (entry->size != 0)
      {
        fprintf(stderr,"existing entry for title \"%s\""
                       "(from file %s, line %d) "
                       "in spacelog entries must have size 0\n",
                       title,filename,line);
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
      gt_spacelog_updateaddentry(entry,filename,line,title,size,work);
    } else
    {
      gt_spacelog_addentry(fcsl,filename,line,title,size,work);
    }
    if (work)
    {
      fcsl->workspace += size;
    } else
    {
      fcsl->splitspace += size;
    }
    gt_firstcodes_updatemax(fcsl);
  } else
  {
    entry = gt_spacelog_find(fcsl,title);
    if (entry == NULL)
    {
      fprintf(stderr,"cannot find title \"%s\" (from file %s, line %d) "
                     "in spacelog entries\n",title,filename,line);
      (void) gt_firstcodes_spacelog_showentries(stderr,fcsl);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    if ((entry->work && !work) || (!entry->work && work))
    {
      fprintf(stderr,"for title \"%s\" (from file %s, line %d) "
                     "in spacelog entries: inconsistent work/splitassignment\n",
               title,filename,line);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    if (work)
    {
      if (entry->size > fcsl->workspace)
      {
        gt_firstcodes_subtract_error(title, filename, line, entry->size, true,
                                     fcsl->workspace);
        (void) gt_firstcodes_spacelog_showentries(stderr,fcsl);
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
      fcsl->workspace -= entry->size;
    } else
    {
      if (entry->size > fcsl->splitspace)
      {
        gt_firstcodes_subtract_error(title, filename, line, entry->size, false,
                                     fcsl->splitspace);
        (void) gt_firstcodes_spacelog_showentries(stderr,fcsl);
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
      fcsl->splitspace -= entry->size;
    }
    if (size > 0)
    {
      gt_spacelog_updateaddentry(entry,filename,line,title,size,work);
      if (work)
      {
        fcsl->workspace += size;
      } else
      {
        fcsl->splitspace += size;
      }
      gt_firstcodes_updatemax(fcsl);
      if (size > entry->size)
      {
        add = true;
        size -= entry->size;
      } else
      {
        size = entry->size - size;
      }
    } else
    {
      size = entry->size;
      entry->size = 0;
    }
  }
  logspace = fcsl->workspace+fcsl->splitspace;
  gt_log_log("file %s, line %d: %s %.2f MB for %s to %sspace; "
             "work=%.2f, split=%.2f, all=%.2f MB (%lu)",
             filename,
             line,
             add ? "add" : "delete",
             GT_MEGABYTES(size),
             title,
             work ? "work" : "split",
             GT_MEGABYTES(fcsl->workspace),
             GT_MEGABYTES(fcsl->splitspace),
             GT_MEGABYTES(logspace),
             (unsigned long) logspace);
  if (gt_ma_enabled() && gt_fa_enabled())
  {
    unsigned long realspace = gt_ma_get_space_current() +
                              gt_fa_get_space_current();
    gt_log_log("current space usage %.2f MB (%.2f+%.2f)",
                                GT_MEGABYTES(realspace),
                                GT_MEGABYTES(gt_ma_get_space_current()),
                                GT_MEGABYTES(gt_fa_get_space_current()));
#ifdef SKDEBUG
    if (fcsl->calc_difference)
    {
      double percent_difference;

      if ((unsigned long) logspace > realspace)
      {
        fprintf(stderr,"overestimating logspace\n");
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
      if (realspace >= 1000000UL)
      {
        /*
        printf("realspace=%lu,logspace=%lu\n",realspace,
                                              (unsigned long) logspace);
        */
        percent_difference = 100.0 * (double) (realspace - logspace)/realspace;
        if (gt_double_larger_double(percent_difference,3.0))
        {
          fprintf(stderr,"space difference of %.2f%% is too large\n",
                       percent_difference);
          exit(GT_EXIT_PROGRAMMING_ERROR);
        }
        if (gt_double_smaller_double(fcsl->max_percent_difference,
                                     percent_difference))
        {
          fcsl->max_percent_difference = percent_difference;
        }
      }
    }
#endif
  }
}
