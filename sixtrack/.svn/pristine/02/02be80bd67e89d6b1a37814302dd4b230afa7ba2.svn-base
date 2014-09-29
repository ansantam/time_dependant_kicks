# rm *.log *.aux
clean:
	rm -rf *.aux *.log *~

# rm also *.pdf
sweepall: clean
	rm -rf *.pdf

all: sweepall
	pdflatex createFlowCharts.tex
	./separatePlots.sh
