// Package gt provides Go bindings for libgenometools.
package gt

/*
#cgo pkg-config: cairo pango
#cgo CFLAGS: -I./../src -I./../obj
#cgo LDFLAGS: -L./../lib -lgenometools

#include "genometools.h"
*/
import "C"

func init() {
	C.gt_lib_init()
}
