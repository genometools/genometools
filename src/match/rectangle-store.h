#ifndef RECTANGLE_STORE_H
#define RECTANGLE_STORE_H
#include <stdint.h>

typedef struct GtArrayGtDiagbandseedRectangle GtArrayGtDiagbandseedRectangle;

typedef struct
{
  uint32_t a_start, a_end, b_start, b_end;
} GtDiagbandseedRectangle;

GtArrayGtDiagbandseedRectangle *gt_rectangle_store_new(void);

void gt_rectangle_store_delete(GtArrayGtDiagbandseedRectangle *rectangle_store);

void gt_rectangle_store_add(GtArrayGtDiagbandseedRectangle *rectangle_store,
                            const GtDiagbandseedRectangle *key);

bool gt_rectangle_overlap(
                   const GtArrayGtDiagbandseedRectangle *rectangle_store,
                   const GtDiagbandseedRectangle *key);

void gt_rectangle_store_show(const GtArrayGtDiagbandseedRectangle
                                *rectangle_store);
#endif
