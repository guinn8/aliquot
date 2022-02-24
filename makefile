PDF_DIR=pdf
COMPILER=pdflatex
.PHONY: pdf

all: 502_Proposal BasicPresentation interm_report AliquotParents TechPresentation Understanding_PollPom

502_Proposal:
	cd tex/$@ && $(COMPILER) $@.tex
	mv tex/$@/$@.pdf pdf

BasicPresentation:
	cd tex/$@ && $(COMPILER) $@.tex
	mv tex/$@/$@.pdf pdf

interm_report:
	cd tex/$@ && $(COMPILER) $@.tex
	mv tex/$@/$@.pdf pdf
	
AliquotParents:
	cd tex/$@ && $(COMPILER) $@.tex
	mv tex/$@/$@.pdf pdf

TechPresentation:
	cd tex/$@ && $(COMPILER) $@.tex
	mv tex/$@/$@.pdf pdf

Understanding_PollPom:
	cd tex/$@ && $(COMPILER) $@.tex
	mv tex/$@/$@.pdf pdf