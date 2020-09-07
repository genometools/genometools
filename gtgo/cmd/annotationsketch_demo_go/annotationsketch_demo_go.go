package main

import (
	"flag"
	"fmt"
	"os"

	gt "github.com/genometools/genometools/gtgo"
)

func sketch(gff3file, pngfile, styleFile string, width uint) error {
	fi := gt.FeatureIndexMemoryNew()
	if err := fi.AddGFF3File(gff3file); err != nil {
		return err
	}
	seqID, err := fi.GetFirstSeqID()
	if err != nil {
		return err
	}
	r, err := fi.GetRangeForSeqID(seqID)
	if err != nil {
		return err
	}
	style, err := gt.StyleNew()
	if err != nil {
		return err
	}
	if err := style.LoadFile(styleFile); err != nil {
		return err
	}
	diagram, err := gt.DiagramNew(fi, seqID, *r, style)
	if err != nil {
		return err
	}
	layout, err := gt.LayoutNew(diagram, width, style)
	if err != nil {
		return err
	}
	height, err := layout.GetHeight()
	if err != nil {
		return err
	}
	canvas, err := gt.CanvasCairoFileNew(style, width, height)
	if err != nil {
		return err
	}
	if err := layout.Sketch(canvas.Canvas()); err != nil {
		return err
	}
	if err := canvas.Write(pngfile); err != nil {
		return err
	}
	return nil
}

func fatal(err error) {
	fmt.Fprintf(os.Stderr, "%s: error: %s\n", os.Args[0], err)
	os.Exit(1)
}

func usage() {
	fmt.Fprintf(os.Stderr, "Usage: %s GFF3_file PNG_file\n", os.Args[0])
	os.Exit(2)
}

func main() {
	styleFile := flag.String("style", "gtdata/sketch/default.style", "Default style file")
	width := flag.Uint("width", 800, "Image width")
	flag.Parse()
	if flag.NArg() != 2 {
		usage()
	}
	if err := sketch(flag.Arg(0), flag.Arg(1), *styleFile, *width); err != nil {
		fatal(err)
	}
}
