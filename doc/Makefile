
CONVERTER := $(shell which rst2html 2>/dev/null || which rst2html.py) --stylesheet=page.css
LATEXCONVERTER := $(shell which rst2latex 2>/dev/null || which rst2latex.py)
PDFCONVERTER := $(shell which rst2pdf 2>/dev/null || which rst2pdf.py) -c
PDFCONCAT := gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite

HTMLDOCFILES = $(patsubst %.rst,%.html,$(wildcard *.rst))
PDFDOCFILES = $(patsubst %.rst,%.pdf,$(wildcard *.rst))

all: manual-html api-doc book.pdf

api-doc: ../doxygen.config
	cd .. && doxygen doxygen.config

manual-html: $(HTMLDOCFILES)

api.pdf: api-doc
	cd api/latex/ && $(MAKE) refman.pdf > make.log
	cp api/latex/refman.pdf api.pdf

book.pdf: $(PDFDOCFILES)
	$(PDFCONCAT) -sOutputFile=book.pdf index.pdf contact.pdf general.pdf manual.pdf faq.pdf others.pdf license.pdf

clean: 
	rm -f $(HTMLDOCFILES)
	rm -f $(PDFDOCFILES)
	rm -rf api

%.html: %.rst header.txt footer.txt
	$(CONVERTER) $< > $@

%.tex: %.rst header.txt footer.txt
	$(LATEXCONVERTER) $< > $@
%.pdf: %.rst header.txt footer.txt
	$(PDFCONVERTER) $< -o $@
#pdflatex -interaction=batchmode $<

.PHONY: all api-doc manual-html clean
