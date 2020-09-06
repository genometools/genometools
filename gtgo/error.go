package gt

/*
#include "genometools.h"
*/
import "C"

import (
	"errors"
	"runtime"
)

// Error object.
type Error struct {
	e *C.GtError
}

// ErrorNew returns a new error object.
func ErrorNew() *Error {
	err := &Error{C.gt_error_new()}
	runtime.SetFinalizer(err, (*Error).delete)
	return err
}

// Get the error string stored in e (the error must be set).
func (e *Error) Get() error {
	return errors.New(C.GoString(C.gt_error_get(e.e)))
}

func (e *Error) delete() {
	C.gt_error_delete(e.e)
}
