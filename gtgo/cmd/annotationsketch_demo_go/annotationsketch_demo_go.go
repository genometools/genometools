package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/genometools/genometools/gtgo"
)

func sketch(filename string) error {
	fi := gt.FeatureIndexMemoryNew()
	if err := fi.AddGFF3File(filename); err != nil {
		return err
	}
	seqID, err := fi.GetFirstSeqID()
	if err != nil {
		return err
	}
	fmt.Println(seqID)

	// TODO

	return nil
}

func fatal(err error) {
	fmt.Fprintf(os.Stderr, "%s: error: %s\n", os.Args[0], err)
	os.Exit(1)
}

func usage() {
	fmt.Fprintf(os.Stderr, "usage: %s GFF3_file\n", os.Args[0])
	os.Exit(2)
}

func main() {
	// TODO flags:
	// - sequence region
	// - range
	// - image width

	flag.Parse()
	if flag.NArg() != 1 {
		usage()
	}
	if err := sketch(flag.Arg(0)); err != nil {
		fatal(err)
	}
}
