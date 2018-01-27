pdflatex main_formated.tex REM -interaction=batchmode
bibtex main_formated
REM makeindex main_formated.nlo -s nomencl.ist -o main_formated.nl
makeindex -s nomencl.ist -t main_formated.nlg -o main_formated.nls main_formated.nlo
pdflatex main_formated.tex -interaction=batchmode
pdflatex main_formated.tex -interaction=batchmode
sumatrapdf main_formated.pdf