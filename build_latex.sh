#!/usr/bin/env bash
# Build paper.pdf from paper.md (run from the project root). Needs pandoc +
# pdflatex. Works on both macOS and Linux: the in-place edits use perl rather
# than sed, because BSD sed (macOS) and GNU sed (Linux) disagree on both the
# -i flag and \n in replacements.
set -e

pandoc paper.md -s -o paper.tex \
  --shift-heading-level-by=-1 \
  -V documentclass=article -V classoption=11pt \
  -V geometry:margin=1in -V fontfamily=mathptmx \
  -H latex_preamble.tex

# LaTeX auto-numbers floats, so drop the manual "Figure N." / "Table X." labels
perl -pi -e 's/\\caption\{\\textbf\{(?:Figure|Table)[^}]*\} /\\caption{/g' paper.tex

# start the appendix and relabel its table as "A1" to match the text
perl -pi -e 's/\\section\{Appendix A\./\\appendix\n\\setcounter{table}{0}\n\\renewcommand{\\thetable}{A\\arabic{table}}\n\\section{Appendix A./' paper.tex

pdflatex -interaction=nonstopmode paper.tex >/dev/null
pdflatex -interaction=nonstopmode paper.tex >/dev/null

if command -v pdfinfo >/dev/null 2>&1; then
  echo "Built paper.pdf ($(pdfinfo paper.pdf | awk '/Pages/{print $2}') pages)"
else
  echo "Built paper.pdf"
fi
