pdflatex main -interaction=batchmode
bibtex main
makeindex main.nlo -s nomencl.ist -o main.nl
pdflatex main.tex -interaction=batchmode
pdflatex main.tex -interaction=batchmode
sumatrapdf main.pdf