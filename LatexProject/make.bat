pdflatex main
bibtex main
makeindex main.nlo -s nomencl.ist -o main.nl
pdflatex main.tex
pdflatex main.tex
sumatrapdf main.pdf