pdflatex main --interaction=batchmode
bibtex main --interaction=batchmode
makeindex main.nlo -s nomencl.ist -o main.nls
pdflatex main.tex --interaction=batchmode
pdflatex main.tex --interaction=batchmode
REM main.pdf
