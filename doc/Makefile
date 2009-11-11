
DOCFILES = $(patsubst %.rst,%.html,$(wildcard *.rst))

all: $(DOCFILES) readme

%.html: %.rst
	rst2html $< > $@

readme.html: readme
	rst2html $< > $@

