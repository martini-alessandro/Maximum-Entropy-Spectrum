all: mlgw.pdf

mlgw.pdf: mlgw.tex mlgw.bbl mlgw.aux
	pdflatex -synctex=1 -interaction=nonstopmode "mlgw".tex

mlgw.aux:
	pdflatex -synctex=1 -interaction=nonstopmode "mlgw".tex

mlgw.bbl: mlgw.aux
	bibtex "mlgw".aux

clean:
	rm -f mlgw.bbl mlgw.out mlgwNotes.bib mlgw.aux mlgw.log mlgw.blg mlgw.toc mlgw.synctex.gz mlgw.pdf
	rm -rf scripts/__pycache__
