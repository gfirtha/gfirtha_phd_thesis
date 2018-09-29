pdflatex booklet_en.tex REM -interaction=batchmode
biber booklet_en
pdflatex booklet_en.tex -interaction=batchmode
pdflatex booklet_en.tex -interaction=batchmode
sumatrapdf booklet_en.pdf