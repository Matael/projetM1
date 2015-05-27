#
# Makefile
# matael, 2015-05-22 16:54
#

TEXC=pdflatex
texc_options=
main_source=main
Bibtex=bibtex

all:
	$(TEXC) $(main_source).tex $(texc_options)
	$(Bibtex) $(main_source)
	$(TEXC) $(main_source).tex $(texc_options)
	$(TEXC) $(main_source).tex $(texc_options)



