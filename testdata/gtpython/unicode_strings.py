#!/usr/bin/python
# -*- coding: utf-8 -*-

from gt.core import *
from gt.extended import *
from gt.annotationsketch import *
from gt.annotationsketch.custom_track import CustomTrack
from gt.core.gtrange import Range
import sys


class GeneRegion(object):
    """Represents a gene region"""

    def __init__(self, start, end):
        self.start = start
        self.end = end

class BetterGeneRegion(GeneRegion):
    """Has even more info and inherits from GeneRegion"""
    def __init__(self, start, end, chromosome, strand):
        GeneRegion.__init__(self, start, end)
        self.chromosome = chromosome
        self.strand = strand

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: " + (sys.argv)[0] +
                         " style_file PNG_file\n")
        sys.exit(1)

    nodes = []

    # construct a gene on the forward strand with two exons
    G = BetterGeneRegion(113558120L, 113578742L, 'chr7', '+')
    seqid = G.chromosome


    gene = FeatureNode.create_new(unicode(seqid),
                                  "gene",
                                   #100,
                                   #900,
                                   G.start,
                                   G.end,
                                   "+")
    gene.add_attribute("Name", u"Testcase™ あがぃいぅ".encode("utf-8"))

    nodes.append(gene)

    exon = FeatureNode.create_new(seqid, "exon", 100, 200, "+")
    gene.add_child(exon)
    intron = FeatureNode.create_new(seqid, "intron", 201, 799, "+")
    gene.add_child(intron)
    exon = FeatureNode.create_new(seqid, "exon", 800, 900, "+")
    gene.add_child(exon)

  # construct a single-exon gene on the reverse strand
  # (within the intron of the forward strand gene)

    reverse_gene = FeatureNode.create_new(seqid, "gene", 400, 600, "-")
    reverse_exon = FeatureNode.create_new(seqid, "exon", 400, 600, "-")
    reverse_gene.add_child(reverse_exon)

    nodes.append(reverse_gene)

    pngfile = (sys.argv)[2]

    style = Style()
    style.load_file((sys.argv)[1])

    diagram = Diagram.from_array(nodes, Range(G.start, G.end),
                                 style)

    layout = Layout(diagram, 600, style)
    height = layout.get_height()
    canvas = CanvasCairoFile(style, 600, height)
    layout.sketch(canvas)

    canvas.to_file(pngfile)
