package gt

/*
#include "genometools.h"
*/
import "C"

import (
	"unsafe"
)

func GFF3InStreamNew(filename string, sorted bool) *NodeStream {
	var (
		ns        *C.GtNodeStream
		cFilename *C.char
	)
	if filename != "" {
		cFilename = C.CString(filename)
	}
	if sorted {
		ns = C.gt_gff3_in_stream_new_sorted(cFilename)
	} else {
		if filename != "" {
			ns = C.gt_gff3_in_stream_new_unsorted(1, &cFilename)
		} else {
			ns = C.gt_gff3_in_stream_new_unsorted(0, nil)
		}
	}
	C.free(unsafe.Pointer(cFilename))
	return nodeStreamNew(ns)
}
