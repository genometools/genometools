ALL=ltrharvest.pdf\
    uniquesub.pdf\
    matstat.pdf\
    mgth.pdf\
    packedindex.pdf\
    tallymer.pdf

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

read1.fna:${GTTESTDATA}/Iowa/at1MB
	head -n 6999 $< > $@

read2.fna:${GTTESTDATA}/Iowa/at1MB
	echo ">" > $@
	tail -n 6945 $< >> $@

U89959.fna:../../testdata/U89959_genomic.fas
	cp $< $@

${TALLYMERTMP}.tex:tallymer.tex read1.fna read2.fna U89959.fna
	./findemptyoutfiles.sh
	./execcommand.pl $< > $@

tallymer.pdf:${TALLYMERTMP}.tex
	pdflatex ${TALLYMERTMP}
	pdflatex ${TALLYMERTMP}
	mv ${TALLYMERTMP}.pdf $@

tallymer.dvi:${TALLYMERTMP}.tex
	latex ${TALLYMERTMP}
	latex ${TALLYMERTMP}
	mv ${TALLYMERTMP}.dvi $@

%.pdf:%.tex
	pdflatex $<
	pdflatex $<

%.dvi:%.tex
	latex $<
	latex $<

clean:
	${RM} *.aux  *.log *.bbl  *.blg *.dvi *.mer *.mct *.mbd comment.cut
	${RM} U89959.fna read[12].fna 
	${RM} reads.al1 reads.esq reads.lcp reads.llv reads.prj reads.suf
	${RM} ${TALLYMERTMP}.tex

cleanup:clean
	${RM} ${ALL}
