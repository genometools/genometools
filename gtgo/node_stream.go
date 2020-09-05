package gt

/*
#include "genometools.h"
*/
import "C"

import (
	"runtime"
)

type NodeStream struct {
	ns  *C.GtNodeStream
	err *Error
}

func nodeStreamNew(ns *C.GtNodeStream) *NodeStream {
	r := &NodeStream{ns, ErrorNew()}
	runtime.SetFinalizer(r, (*NodeStream).delete)
	return r
}

func (ns *NodeStream) Next() (*GenomeNode, error) {
	var gn *C.GtGenomeNode
	if C.gt_node_stream_next(ns.ns, &gn, ns.err.e) != 0 {
		return nil, ns.err.Get()
	}
	if gn == nil {
		return nil, nil
	}
	return genomeNodeNew(gn), nil
}

func (ns *NodeStream) delete() {
	C.gt_node_stream_delete(ns.ns)
}
