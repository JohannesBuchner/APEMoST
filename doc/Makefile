all: manual-html api-doc api.pdf book.pdf

api-doc: ../doxygen.config
	cd .. && doxygen doxygen.config

_build/latex/APEMoST.pdf:
	$(MAKE) -f Makefile.sphinx latexpdf

manual-html: 
	$(MAKE) -f Makefile.sphinx html

api/latex/refman.pdf: api-doc
	cd api/latex/ && $(MAKE) refman.pdf > make.log

website: _build/latex/APEMoST.pdf manual-html api/latex/refman.pdf
	mkdir -p _build/website
	cp -uv _build/latex/APEMoST.pdf _build/website/book.pdf
	cp -uv api/latex/refman.pdf _build/website/api.pdf
	cp -uv -r _build/html/* _build/website/
	mkdir -p _build/website/api/
	cp -uv -r api/html/* _build/website/api/

clean: 
	$(MAKE) -f Makefile.sphinx clean
	rm -f *.html
	rm -f *.pdf
	rm -rf _build/website

.PHONY: all api-doc manual-html clean
