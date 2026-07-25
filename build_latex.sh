#!/usr/bin/env bash
# Build paper.pdf from paper.md (run from the project root). Needs pandoc + pdflatex.
set -e
pandoc paper.md -s -o paper.tex \
  --shift-heading-level-by=-1 \
  -V documentclass=article -V classoption=11pt \
  -V geometry:margin=1in -V fontfamily=mathptmx \
  -H latex_preamble.tex
# LaTeX auto-numbers floats, so drop the manual "Figure N." / "Table X." labels
sed -i -E 's/\\caption\{\\textbf\{(Figure|Table)[^}]*\} /\\caption{/g' paper.tex
# start the appendix and relabel its table as "A1" to match the text
sed -i 's/\\section{Appendix A\./\\appendix\n\\setcounter{table}{0}\n\\renewcommand{\\thetable}{A\\arabic{table}}\n\\section{Appendix A./' paper.tex
pdflatex -interaction=nonstopmode paper.tex >/dev/null
pdflatex -interaction=nonstopmode paper.tex >/dev/null
echo "Built paper.pdf ($(pdfinfo paper.pdf | awk '/Pages/{print $2}') pages)"
