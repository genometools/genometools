ALL=ltrharvestman.pdf\
    uniquesub.pdf\
    matstat.pdf\
    mgth.pdf\
    packedindex.pdf

all:${ALL}

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
	-diff $@ $*.bbl || mv $*.bbl $@ && rm $*.bbl
	mv $*.aux $*-pdf.aux

%-pdf.aux: %.tex
	-test -e $@ && mv $@ $*.aux
	pdflatex $< && rm $*.pdf
	mv $*.aux $@


clean:
	${RM} *.aux  *.log *.bbl  *.blg *.dvi comment.cut

cleanup:clean
	${RM} ${ALL}
