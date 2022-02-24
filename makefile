PDF_DIR=pdf
COMPILER=pdflatex
.PHONY: pdf

all: experiments_densities_kparent_proposal kparent_density_basic_presentation kparent_aliquot_interm_report pollack_pomerance_kparent kparent_density_technical_presentation pollack_pomerance_model_notes

experiments_densities_kparent_proposal:
	cd tex/$@ && $(COMPILER) $@.tex
	mv tex/$@/$@.pdf pdf

kparent_density_basic_presentation:
	cd tex/$@ && $(COMPILER) $@.tex
	mv tex/$@/$@.pdf pdf

kparent_aliquot_interm_report:
	cd tex/$@ && $(COMPILER) $@.tex
	mv tex/$@/$@.pdf pdf
	
pollack_pomerance_kparent:
	cd tex/$@ && $(COMPILER) $@.tex
	mv tex/$@/$@.pdf pdf

kparent_density_technical_presentation:
	cd tex/$@ && $(COMPILER) $@.tex
	mv tex/$@/$@.pdf pdf

pollack_pomerance_model_notes:
	cd tex/$@ && $(COMPILER) $@.tex
	mv tex/$@/$@.pdf pdf