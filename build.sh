#!/bin/bash
pdflatex grant_report.tex
bibtex grant_report.aux
pdflatex grant_report.tex
pdflatex grant_report.tex
