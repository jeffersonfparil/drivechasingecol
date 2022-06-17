# Build manuscript
```{sh}
cd ms/
pdflatex main.tex --bibliography references.bib -o main.pdf
bibtex main
pdflatex main.tex --bibliography references.bib -o main.pdf
pdflatex main.tex --bibliography references.bib -o main.pdf
```

# Track changes between revisions
```{sh}
latexdiff \
    main.tex \
    ~/Documents/drivechasingecol/ms/main.tex \
    > tracked_changes_with_latexdiff_20220617.tex
```