package gt

/*
#include "genometools.h"
*/
import "C"

import (
	"errors"
	"runtime"
)

type Error struct {
	e *C.GtError
}

func ErrorNew() *Error {
	err := &Error{C.gt_error_new()}
	runtime.SetFinalizer(err, (*Error).delete)
	return err
}

func (e *Error) Get() error {
	return errors.New(C.GoString(C.gt_error_get(e.e)))
}

func (e *Error) delete() {
	C.gt_error_delete(e.e)
}
