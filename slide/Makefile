# Makefile para compilar apresentação Beamer

MAIN = apresentacao
TEX = $(MAIN).tex
PDF = $(MAIN).pdf

# Compilador
LATEX = pdflatex
BIBTEX = bibtex

# Flags
LATEX_FLAGS = -interaction=nonstopmode -halt-on-error

.PHONY: all clean view cleanall

all: $(PDF)

$(PDF): $(TEX)
	$(LATEX) $(LATEX_FLAGS) $(TEX)
	$(LATEX) $(LATEX_FLAGS) $(TEX)

clean:
	rm -f *.aux *.log *.nav *.out *.snm *.toc *.vrb

cleanall: clean
	rm -f $(PDF)

view: $(PDF)
	open $(PDF) || xdg-open $(PDF) || evince $(PDF)

# Compilação rápida (1 passagem)
quick:
	$(LATEX) $(LATEX_FLAGS) $(TEX)

# Compilação completa (2 passagens para sumário)
full: all
