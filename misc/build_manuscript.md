# Track changes between revisions
```{sh}
latexdiff \
    ~/Documents/drivechasingecol/ms/main.tex \
    main.tex \
    > tracked_changes_with_latexdiff_20220617.tex
```

# Build manuscript
```{sh}
pdflatex -interaction=nonstopmode main.tex  --bibliography references.bib -o main.pdf
bibtex main
pdflatex -interaction=nonstopmode main.tex  --bibliography references.bib -o main.pdf
pdflatex -interaction=nonstopmode main.tex  --bibliography references.bib -o main.pdf

pdflatex -interaction=nonstopmode tracked_changes_with_latexdiff_20220617.tex  --bibliography references.bib -o tracked_changes_with_latexdiff_20220617.pdf
bibtex tracked_changes_with_latexdiff_20220617
pdflatex -interaction=nonstopmode tracked_changes_with_latexdiff_20220617.tex  --bibliography references.bib -o tracked_changes_with_latexdiff_20220617.pdf
pdflatex -interaction=nonstopmode tracked_changes_with_latexdiff_20220617.tex  --bibliography references.bib -o tracked_changes_with_latexdiff_20220617.pdf

pdflatex -interaction=nonstopmode Responses_to_reviewers_20220617.tex  --bibliography references.bib -o Responses_to_reviewers_20220617.pdf
bibtex Responses_to_reviewers_20220617
pdflatex -interaction=nonstopmode Responses_to_reviewers_20220617.tex  --bibliography references.bib -o Responses_to_reviewers_20220617.pdf
pdflatex -interaction=nonstopmode Responses_to_reviewers_20220617.tex  --bibliography references.bib -o Responses_to_reviewers_20220617.pdf

rm *.aux *.bbl *.blg *.fdb_latexmk *.fls *.log *.out *.xml *-blx.bib
```
