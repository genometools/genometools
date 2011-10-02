#include <string.h>
#include "core/ma.h"
#include "core/fa.h"
#include "core/log_api.h"
#include "core/spacecalc.h"
#include "core/spacepeak.h"
#include "firstcodes-spacelog.h"

typedef struct
{
  const char *filename;
  int line;
  const char *title;
  size_t size;
} GtFirstcodespacelogentry;

struct GtFirstcodesspacelog
{
  size_t workspace,
         splitspace;
  GtFirstcodespacelogentry *entries;
  unsigned long nextfree, allocated;
};

GtFirstcodesspacelog *gt_firstcodes_spacelog_new(void)
{
  GtFirstcodesspacelog *fcsl;

  fcsl = gt_malloc(sizeof (GtFirstcodesspacelog));
  fcsl->workspace = fcsl->splitspace = 0;
  fcsl->nextfree = fcsl->allocated = 0;
  fcsl->entries = NULL;
  return fcsl;
}

size_t gt_firstcodes_spacelog_total(GtFirstcodesspacelog *fcsl)
{
  return fcsl->workspace + fcsl->splitspace;
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
                                       const char *title,size_t size)
{
  entry->filename = filename;
  entry->line = line;
  entry->title = title;
  entry->size = size;
}

static void gt_spacelog_addentry(GtFirstcodesspacelog *fcsl,
                                 const char *filename,int line,
                                 const char *title,size_t size)
{
  gt_assert(fcsl->nextfree <= fcsl->allocated);
  if (fcsl->nextfree == fcsl->allocated)
  {
    fcsl->allocated += 16UL;
    fcsl->entries = gt_realloc(fcsl->entries,sizeof (*fcsl->entries) *
                                             fcsl->allocated);
  }
  gt_spacelog_updateaddentry(fcsl->entries + fcsl->nextfree,
                             filename,line,title,size);
  fcsl->nextfree++;
}

static bool gt_spacelog_showentries(FILE *fp,const GtFirstcodesspacelog *fcsl)
{
  unsigned long idx;
  bool foundnonempty = false;

  for (idx = 0; idx < fcsl->nextfree; idx++)
  {
    if (fcsl->entries[idx].size > 0)
    {
      fprintf(fp,"%s %d %s %lu\n",
               fcsl->entries[idx].filename,
               fcsl->entries[idx].line,
               fcsl->entries[idx].title,
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
    if (gt_spacelog_showentries(stderr,fcsl))
    {
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
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

void gt_firstcodes_spacelog_add(GtFirstcodesspacelog *fcsl,
                                int line,
                                const char *filename,
                                bool add,
                                const char *title,
                                bool work,
                                size_t size)
{
  GtFirstcodespacelogentry *entry;

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
      gt_spacelog_updateaddentry(entry,filename,line,title,size);
    } else
    {
      gt_spacelog_addentry(fcsl,filename,line,title,size);
    }
    if (work)
    {
      fcsl->workspace += size;
    } else
    {
      fcsl->splitspace += size;
    }
  } else
  {
    gt_assert(size == 0);
    entry = gt_spacelog_find(fcsl,title);
    if (entry == NULL)
    {
      fprintf(stderr,"cannot find title \"%s\" (from file %s, line %d) "
                     "in spacelog entries\n",title,filename,line);
      (void) gt_spacelog_showentries(stderr,fcsl);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    size = entry->size;
    entry->size = 0;
    if (work)
    {
      if (size > fcsl->workspace)
      {
        gt_firstcodes_subtract_error(title, filename, line, size, true,
                                     fcsl->workspace);
        (void) gt_spacelog_showentries(stderr,fcsl);
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
      fcsl->workspace -= size;
    } else
    {
      if (size > fcsl->splitspace)
      {
        gt_firstcodes_subtract_error(title, filename, line, size, false,
                                     fcsl->splitspace);
        (void) gt_spacelog_showentries(stderr,fcsl);
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
      fcsl->splitspace -= size;
    }
  }
  gt_log_log("file %s, line %d: %s %.2f MB for %s to %sspace; work=%.2f, "
             "split=%.2f, all=%.2f MB",
             filename,
             line,
             add ? "add" : "delete",
             GT_MEGABYTES(size),
             title,
             work ? "work" : "split",
             GT_MEGABYTES(fcsl->workspace),
             GT_MEGABYTES(fcsl->splitspace),
             GT_MEGABYTES(fcsl->workspace+fcsl->splitspace));
  if (gt_ma_enabled() && gt_fa_enabled())
  {
    gt_log_log("current space usage %.2f",
               GT_MEGABYTES(gt_ma_get_space_current() +
                            gt_fa_get_space_current()));
  }
}
