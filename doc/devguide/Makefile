latexfile = devguide

TEX = pdflatex

$(latexfile).pdf : $(latexfile).tex
	while ($(TEX) --shell-escape $(latexfile) ; \
	grep -q "Rerun to get cross" $(latexfile).log ) do true ; \
	done


pdf : $(latexfile).pdf

$(latexfile).ps : $(latexfile).dvi
	dvips -o  $(latexfile).ps $(latexfile).dvi

ps : $(latexfile).ps 

clean: 
	rm -f $(latexfile).pdf *.aux *.ps *.dvi *.log *.out *.toc *.bak
