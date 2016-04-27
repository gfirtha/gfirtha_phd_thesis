pdflatex main
bibtex main
makeindex main.nlo -s nomencl.ist -o main.nls
pdflatex main.tex
main.pdf