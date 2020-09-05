package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/genometools/genometools/gtgo"
)

func validate(filename string) error {
	// create GFF3 input stream
	inStream := gt.GFF3InStreamNew(filename, false)

	// pull the features through the stream
	for {
		gn, err := inStream.Next()
		if err != nil {
			return err
		}
		if gn == nil {
			return nil
		}
	}
}

func fatal(err error) {
	fmt.Fprintf(os.Stderr, "%s: error: %s\n", os.Args[0], err)
	os.Exit(1)
}

func usage() {
	fmt.Fprintf(os.Stderr, "usage: %s [GFF3_file]\n", os.Args[0])
	os.Exit(2)
}

func main() {
	flag.Parse()
	if flag.NArg() > 1 {
		usage()
	}
	var filename string
	if flag.NArg() == 1 {
		filename = flag.Arg(0)
	}
	if err := validate(filename); err != nil {
		fatal(err)
	}
	fmt.Println("input is valid GFF3")
}
