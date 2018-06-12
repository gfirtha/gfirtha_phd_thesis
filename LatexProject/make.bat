pdflatex main.tex REM -interaction=batchmode
bibtex main
REM makeindex main.nlo -s nomencl.ist -o main.nl
makeindex -s nomencl.ist -t main.nlg -o main.nls main.nlo
pdflatex main_.tex -interaction=batchmode
pdflatex main.tex -interaction=batchmode
sumatrapdf main.pdf