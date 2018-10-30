pdflatex main.tex
bibtex main
REM makeindex main.nlo -s nomencl.ist -o main.nl
makeindex -s nomencl.ist -t main.nlg -o main.nls main.nlo
pdflatex main.tex -interaction=batchmode
pdflatex main.tex 2>error.txt >output.txt
sumatrapdf main.pdf