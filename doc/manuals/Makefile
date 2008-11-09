ALL=ltrharvest.pdf\
    uniquesub.pdf\
    matstat.pdf\
    mgth.pdf\
    packedindex.pdf

TALLYMERTMP=tallymer-tmp

all:${ALL}

annotationsketch: annotationsketch.pdf

%.pdf:%.tex %-pdf.bbl
	-test -e $*-pdf.aux && mv $*-pdf.aux $*.aux
	-test -e $*-pdf.bbl && mv $*-pdf.bbl $*.bbl
	pdflatex $<
	pdflatex $<
	mv $*.aux $*-pdf.aux
	-test -e $*.bbl && mv $*.bbl $*-pdf.bbl

.PRECIOUS: %-pdf.bbl %-dvi.bbl %-dvi.aux %-pdf.aux
%-pdf.bbl: gtmanuals.bib %-pdf.aux
	-test -e $*-pdf.aux && mv $*-pdf.aux $*.aux
	-bibtex $*
	test -e $@ || touch $@
	-diff $@ $*.bbl || mv $*.bbl $@ && rm -f $*.bbl
	mv $*.aux $*-pdf.aux

%-pdf.aux: %.tex
	-test -e $@ && mv $@ $*.aux
	pdflatex $< && rm $*.pdf
	mv $*.aux $@

${TALLYMERTMP}.tex:tallymer.tex
	./findemptyoutfiles.sh
	head -n 6999 ${GTTESTDATA}/Iowa/at1MB > read1.fna
	echo ">" > read2.fna
	tail -n 6946 ${GTTESTDATA}/Iowa/at1MB > read2.fna
	cp ../../testdata/U89959_genomic.fas U89959.fna
	./execcommand.pl $< > $@

tallymer.pdf:${TALLYMERTMP}.tex
	pdflatex ${TALLYMERTMP}
	pdflatex ${TALLYMERTMP}
	mv ${TALLYMERTMP}.pdf $@

tallymer.dvi:${TALLYMERTMP}.tex
	latex ${TALLYMERTMP}
	latex ${TALLYMERTMP}
	mv ${TALLYMERTMP}.dvi $@

clean:
	${RM} *.aux  *.log *.bbl  *.blg *.dvi *.mer *.mct *.mbd comment.cut
	${RM} U89959.fna read[12].fna 
	${RM} reads.al1 reads.esq reads.lcp reads.llv reads.prj reads.suf
	${RM} ${TALLYMERTMP}.tex

cleanup:clean
	${RM} ${ALL}
