# Build manuscript
```{sh}
cd ms/
pdflatex main.tex --bibliography references.bib -o main.pdf
bibtex main
pdflatex main.tex --bibliography references.bib -o main.pdf
pdflatex main.tex --bibliography references.bib -o main.pdf
```
