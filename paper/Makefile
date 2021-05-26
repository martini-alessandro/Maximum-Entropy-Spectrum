all: mesa.pdf

mesa.pdf: mesa.tex mesa.bbl mesa.aux Bibliography.bib
	pdflatex -synctex=1 -interaction=nonstopmode "mesa".tex

mesa.aux: Bibliography.bib
	pdflatex -synctex=1 -interaction=nonstopmode "mesa".tex

mesa.bbl: mesa.aux  Bibliography.bib
	bibtex "mesa".aux

clean:
	rm -f mesa.bbl mesa.out mesaNotes.bib mesa.aux mesa.log mesa.blg mesa.toc mesa.synctex.gz mesa.pdf
	rm -rf scripts/__pycache__
