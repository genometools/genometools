// +build cairo

package gt

/*
#include "genometools.h"
*/
import "C"

import (
	"runtime"
	"unsafe"
)

// Diagram object.
type Diagram struct {
	d *C.GtDiagram
}

// DiagramNew returns a new diagram object representing the feature nodes in fi
// in region seqID overlapping with r. The style object will be used to
// determine collapsing options during the layout process.
func DiagramNew(fi *FeatureIndex, seqID string, r Range, style *Style) (*Diagram, error) {
	e := ErrorNew()
	s := C.CString(seqID)
	defer C.free(unsafe.Pointer(s))
	var cr C.GtRange
	cr.start = C.ulong(r.Start)
	cr.end = C.ulong(r.End)
	res := C.gt_diagram_new(fi.fi, s, &cr, style.s, e.e)
	if res == nil {
		return nil, e.Get()
	}
	d := &Diagram{res}
	runtime.SetFinalizer(d, (*Diagram).delete)
	return d, nil
}

func (d *Diagram) delete() {
	C.gt_diagram_delete(d.d)
}
