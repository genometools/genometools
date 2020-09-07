package gt

/*
#include "genometools.h"
*/
import "C"

import (
	"runtime"
	"unsafe"
)

// FeatureIndex interface.
type FeatureIndex struct {
	fi *C.GtFeatureIndex
}

// FeatureIndexMemoryNew creates a new memory feature index object.
func FeatureIndexMemoryNew() *FeatureIndex {
	r := &FeatureIndex{C.gt_feature_index_memory_new()}
	runtime.SetFinalizer(r, (*FeatureIndex).delete)
	return r
}

// AddGFF3File adds all features contained in gff3file to fi, if gff3file is
// valid. Otherwise, fi is not changed and an error is returned.
func (fi *FeatureIndex) AddGFF3File(gff3file string) error {
	e := ErrorNew()
	fn := C.CString(gff3file)
	defer C.free(unsafe.Pointer(fn))
	if C.gt_feature_index_add_gff3file(fi.fi, fn, e.e) != 0 {
		return e.Get()
	}
	return nil
}

// GetFirstSeqID returns the first sequence region identifier added to fi.
func (fi *FeatureIndex) GetFirstSeqID() (string, error) {
	e := ErrorNew()
	s := C.gt_feature_index_get_first_seqid(fi.fi, e.e)
	if s == nil {
		return "", e.Get()
	}
	defer C.free(unsafe.Pointer(s))
	return C.GoString(s), nil
}

// GetSeqIDs returns all sequence region identifiers contained in fi (in
// alphabetical order).
func (fi *FeatureIndex) GetSeqIDs() ([]string, error) {
	e := ErrorNew()
	r := C.gt_feature_index_get_seqids(fi.fi, e.e)
	if r == nil {
		return nil, e.Get()
	}
	sa := strArrayNew(r)
	return sa.GetArray(), nil
}

// GetRangeForSeqID returns the range of all features contained in fi for the
// region identifier seqid.
func (fi *FeatureIndex) GetRangeForSeqID(seqID string) (*Range, error) {
	var r C.GtRange
	s := C.CString(seqID)
	defer C.free(unsafe.Pointer(s))
	e := ErrorNew()
	if C.gt_feature_index_get_range_for_seqid(fi.fi, &r, s, e.e) != 0 {
		return nil, e.Get()
	}
	return &Range{uint(r.start), uint(r.end)}, nil
}

func (fi *FeatureIndex) delete() {
	C.gt_feature_index_delete(fi.fi)
}
