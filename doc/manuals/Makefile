ALL=ltrharvestman.pdf\
    uniquesub.pdf

all:${ALL}

%.pdf:%.tex
	pdflatex $<
	pdflatex $<

clean:
	${RM} *.aux  *.log *.bbl  *.blg *.dvi comment.cut

cleanup:clean
	${RM} ${ALL}
