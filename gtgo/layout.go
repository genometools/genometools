// +build cairo

package gt

/*
#include "genometools.h"
*/
import "C"

import (
	"runtime"
)

// Layout object.
type Layout struct {
	l *C.GtLayout
}

// LayoutNew creates a new Layout object for the contents of d. The layout is
// done for a target image of given width and using the rules in s.
func LayoutNew(d *Diagram, width uint, s *Style) (*Layout, error) {
	e := ErrorNew()
	r := C.gt_layout_new(d.d, C.uint(width), s.s, e.e)
	if r == nil {
		return nil, e.Get()
	}
	l := &Layout{r}
	runtime.SetFinalizer(l, (*Layout).delete)
	return l, nil
}

// GetHeight returns the height of l in pixels.
func (l *Layout) GetHeight() (uint, error) {
	e := ErrorNew()
	var result C.ulong
	if C.gt_layout_get_height(l.l, &result, e.e) != 0 {
		return 0, e.Get()
	}
	return uint(result), nil
}

// Sketch renders c on the target canvas c.
func (l *Layout) Sketch(c *Canvas) error {
	e := ErrorNew()
	if C.gt_layout_sketch(l.l, c.c, e.e) != 0 {
		return e.Get()
	}
	return nil
}

func (l *Layout) delete() {
	C.gt_layout_delete(l.l)
}
