package main

import (
	"flag"
	"fmt"
	"os"

	gt "github.com/genometools/genometools/gtgo"
)

func sketch(filename, styleFile string) error {
	fi := gt.FeatureIndexMemoryNew()
	if err := fi.AddGFF3File(filename); err != nil {
		return err
	}
	seqID, err := fi.GetFirstSeqID()
	if err != nil {
		return err
	}
	fmt.Println(seqID)
	r, err := fi.GetRangeForSeqID(seqID)
	if err != nil {
		return err
	}
	fmt.Printf("%d, %d\n", r.Start, r.End)

	seqIDs, err := fi.GetSeqIDs()
	if err != nil {
		return err
	}
	for _, seqID = range seqIDs {
		fmt.Println(seqID)
	}
	style, err := gt.StyleNew()
	if err != nil {
		return err
	}
	if err := style.LoadFile(styleFile); err != nil {
		return err
	}

	// TODO

	return nil
}

func fatal(err error) {
	fmt.Fprintf(os.Stderr, "%s: error: %s\n", os.Args[0], err)
	os.Exit(1)
}

func usage() {
	fmt.Fprintf(os.Stderr, "Usage: %s GFF3_file\n", os.Args[0])
	os.Exit(2)
}

func main() {
	// TODO flags:
	// - sequence region
	// - range
	// - image width
	styleFile := flag.String("style", "gtdata/sketch/default.style", "Default style file")
	flag.Parse()
	if flag.NArg() != 1 {
		usage()
	}
	if err := sketch(flag.Arg(0), *styleFile); err != nil {
		fatal(err)
	}
}
