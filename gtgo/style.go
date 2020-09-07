package gt

/*
#include "genometools.h"
*/
import "C"

import (
	"runtime"
	"unsafe"
)

// Style object.
type Style struct {
	s *C.GtStyle
}

// StyleNew creates a new style object.
func StyleNew() (*Style, error) {
	e := ErrorNew()
	r, err := C.gt_style_new(e.e)
	if err != nil {
		return nil, err
	}
	s := &Style{r}
	runtime.SetFinalizer(s, (*Style).delete)
	return s, nil
}

// LoadFile loads and executes Lua style file with given filename.
// This file must define a global table called style.
func (s *Style) LoadFile(filename string) error {
	e := ErrorNew()
	fn := C.CString(filename)
	defer C.free(unsafe.Pointer(fn))
	if C.gt_style_load_file(s.s, fn, e.e) != 0 {
		return e.Get()
	}
	return nil
}

func (s *Style) delete() {
	C.gt_style_delete(s.s)
}
