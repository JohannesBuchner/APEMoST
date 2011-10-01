
DOCFILES = $(patsubst %.rst,%.html,$(wildcard *.rst))
CONVERTER = $(shell which rst2html 2>/dev/null || which rst2html.py)


all: $(DOCFILES)

clean: 
	rm $(DOCFILES)

%.html: %.rst header.txt footer.txt
	$(CONVERTER) $< > $@

readme.html: readme
	$(CONVERTER) $< > $@

