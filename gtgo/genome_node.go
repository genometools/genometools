package gt

/*
#include "genometools.h"
*/
import "C"

import (
	"runtime"
)

type GenomeNode struct {
	gn *C.GtGenomeNode
}

func genomeNodeNew(gn *C.GtGenomeNode) *GenomeNode {
	r := &GenomeNode{gn}
	runtime.SetFinalizer(r, (*GenomeNode).delete)
	return r
}

func (gn *GenomeNode) delete() {
	C.gt_genome_node_delete(gn.gn)
}
