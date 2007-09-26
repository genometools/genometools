BASE=ltrharvestman

all:${BASE}.pdf

%.pdf:%.tex
	pdflatex $<
	pdflatex $<

clean:
	${RM} ${BASE}.aux  ${BASE}.log ${BASE}.bbl  ${BASE}.blg ${BASE}.dvi
