# Track changes between revisions
```{sh}
cd ~/Documents/drivechasingecol/ms/
latexdiff \
    main_20220510.tex \
    main.tex \
    > tracked_changes_with_latexdiff.tex
```

# Build manuscript
```{sh}
cd ~/Documents/drivechasingecol/ms/
pdflatex -interaction=nonstopmode main.tex  --bibliography references.bib -o main.pdf
bibtex main
pdflatex -interaction=nonstopmode main.tex  --bibliography references.bib -o main.pdf
pdflatex -interaction=nonstopmode main.tex  --bibliography references.bib -o main.pdf

pdflatex -interaction=nonstopmode tracked_changes_with_latexdiff.tex  --bibliography references.bib -o tracked_changes_with_latexdiff.pdf
bibtex tracked_changes_with_latexdiff
pdflatex -interaction=nonstopmode tracked_changes_with_latexdiff.tex  --bibliography references.bib -o tracked_changes_with_latexdiff.pdf
pdflatex -interaction=nonstopmode tracked_changes_with_latexdiff.tex  --bibliography references.bib -o tracked_changes_with_latexdiff.pdf

pdflatex -interaction=nonstopmode Responses_to_reviewers_round2.tex  --bibliography references.bib -o Responses_to_reviewers_round2.pdf
bibtex Responses_to_reviewers_round2
pdflatex -interaction=nonstopmode Responses_to_reviewers_round2.tex  --bibliography references.bib -o Responses_to_reviewers_round2.pdf
pdflatex -interaction=nonstopmode Responses_to_reviewers_round2.tex  --bibliography references.bib -o Responses_to_reviewers_round2.pdf

rm *.aux *.bbl *.blg *.fdb_latexmk *.fls *.log *.out *.xml *-blx.bib
```

# Count the number of words
```{sh}
cd ~/Documents/drivechasingecol/ms/
detex main.tex | wc -w
```