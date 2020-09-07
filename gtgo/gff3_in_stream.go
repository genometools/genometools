package gt

/*
#include "genometools.h"
*/
import "C"

import (
	"unsafe"
)

// GFF3InStreamNew returns a new GFF3 input stream which subsequently reads
// the GFF3 file denoted by filename.
// If sorted is true, the GFF3 file has to be sorted.
// If the GFF3 file is sorted the memory footprint is O(1) on average.
// Otherwise the memory footprint is O(file size) in the worst-case.
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
