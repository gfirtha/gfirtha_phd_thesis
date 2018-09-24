pdflatex main.tex
bibtex main
REM makeindex main.nlo -s nomencl.ist -o main.nl
makeindex -s nomencl.ist -t main.nlg -o main.nls main.nlo
pdflatex main_.tex
pdflatex main.tex
sumatrapdf main.pdf