// +build cairo

package gt

/*
#include "genometools.h"

static GtCanvasCairoFile* to_canvas_cairo_file_ptr(GtCanvas *c)
{
  return (GtCanvasCairoFile*) c;
}
*/
import "C"

import (
	"unsafe"
)

// CanvasCairoFile object.
// Implements the Canvas interface.
type CanvasCairoFile struct {
	c *Canvas
}

// CanvasCairoFileNew creates a new CanvasCairOFile object.
func CanvasCairoFileNew(s *Style, width, height uint) (*CanvasCairoFile, error) {
	e := ErrorNew()
	r := C.gt_canvas_cairo_file_new(s.s, C.GT_GRAPHICS_PNG, C.ulong(width),
		C.ulong(height), nil, e.e)
	if r == nil {
		return nil, e.Get()
	}
	return &CanvasCairoFile{canvasNew(r)}, nil
}

// Write renders s to the file with filename.
func (c *CanvasCairoFile) Write(filename string) error {
	e := ErrorNew()
	fn := C.CString(filename)
	defer C.free(unsafe.Pointer(fn))
	if C.gt_canvas_cairo_file_to_file(C.to_canvas_cairo_file_ptr(c.c.c), fn, e.e) != 0 {
		return e.Get()
	}
	return nil
}

// Canvas returns the underlying canvas.
func (c *CanvasCairoFile) Canvas() *Canvas {
	return c.c
}
