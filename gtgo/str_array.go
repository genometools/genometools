package gt

/*
#include "genometools.h"
*/
import "C"

import (
	"runtime"
)

// StrArray object.
type StrArray struct {
	sa *C.GtStrArray
}

func strArrayNew(sa *C.GtStrArray) *StrArray {
	r := &StrArray{sa}
	runtime.SetFinalizer(r, (*StrArray).delete)
	return r
}

// Get returns string with number strnum from sa.
// strnum must be smaller than sa.Size().
func (sa *StrArray) Get(strnum uint) string {
	return C.GoString(C.gt_str_array_get(sa.sa, C.ulong(strnum)))
}

// Size returns the number of strings stored in sa.
func (sa *StrArray) Size() uint {
	return uint(C.gt_str_array_size(sa.sa))
}

// GetArray returns a Go string array of all strings in sa.
func (sa *StrArray) GetArray() []string {
	var s []string
	size := sa.Size()
	for i := uint(0); i < size; i++ {
		s = append(s, sa.Get(i))
	}
	return s
}

func (sa *StrArray) delete() {
	C.gt_str_array_delete(sa.sa)
}
