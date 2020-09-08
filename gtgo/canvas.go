// +build cairo

package gt

/*
#include "genometools.h"
*/
import "C"

import (
	"runtime"
)

// Canvas interface.
type Canvas struct {
	c *C.GtCanvas
}

func canvasNew(c *C.GtCanvas) *Canvas {
	r := &Canvas{c}
	runtime.SetFinalizer(r, (*Canvas).delete)
	return r
}

func (c *Canvas) delete() {
	C.gt_canvas_delete(c.c)
}
